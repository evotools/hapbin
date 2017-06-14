/*
 * MPIRPC: MPI based invocation of functions on other ranks
 * Copyright (C) 2014-2017 Colin MacLean <cmaclean@illinois.edu>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "manager.hpp"
#include "common.hpp"
#include <mpi.h>

#define BUFFER_SIZE 10*1024*1024

namespace mpirpc
{

Manager::Manager(MPI_Comm comm) : m_comm(comm), m_nextTypeId(0), m_count(0), m_shutdown(false)
{
    MPI_Comm_rank(m_comm, &m_rank);
    MPI_Comm_size(comm, &m_numProcs);

    const int nitems = 2;
    int blocklengths[2] = {1,1};
    MPI_Datatype types[2] = {MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG};
    MPI_Aint offsets[2];
    offsets[0] = offsetof(ObjectInfo, type);
    offsets[1] = offsetof(ObjectInfo, id);
    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &MpiObjectInfo);
    MPI_Type_commit(&MpiObjectInfo);
    void *buffer = malloc(BUFFER_SIZE);
    MPI_Buffer_attach(buffer, BUFFER_SIZE);
    MPI_Barrier(m_comm);
}

Manager::~Manager()
{
    MPI_Type_free(&MpiObjectInfo);
    for (auto i : m_mpiMessages)
        delete i.second;
    for (auto i : m_mpiObjectMessages)
        i.second.reset();
    for (auto i : m_registeredFunctions)
        delete i.second;
    for (auto i : m_registeredObjects)
        delete i;
}

int Manager::rank() const
{
    return m_rank;
}

void Manager::notifyNewObject(TypeId type, ObjectId id)
{
    if (m_shutdown)
        return;
    std::shared_ptr<ObjectInfo> info(new ObjectInfo(type, id));
    for(int i = 0; i < m_numProcs; ++i)
    {
        if (i != m_rank)
        {
            MPI_Request req;
            MPI_Issend(info.get(), 1, MpiObjectInfo, i, MPIRPC_TAG_NEW, m_comm, &req);
            m_mpiObjectMessages[req] = info;
        }
    }
}

void Manager::sendRawMessage(int rank, const std::vector<char> *data, int tag)
{
    if (checkSends() && !m_shutdown) {
        MPI_Request req;
        MPI_Issend((void*) data->data(), data->size(), MPI_CHAR, rank, tag, m_comm, &req);
        m_mpiMessages[req] = data;
    }
}

void Manager::sendRawMessageToAll(const std::vector<char>* data, int tag)
{
    for (int i = 0; i < m_numProcs; ++i) {
#ifndef USE_MPI_LOCALLY
        if (i != m_rank) {
#endif
            sendRawMessage(i, data, tag);
#ifndef USE_MPI_LOCALLY
        }
#endif
    }
}

void Manager::registerUserMessageHandler(int tag, UserMessageHandler callback) {
    m_userMessageHandlers[tag] = callback;
}

bool Manager::checkSends() {
    for (auto i = m_mpiObjectMessages.begin(); i != m_mpiObjectMessages.end();) {
        MPI_Request req = i->first;
        int flag = true;
        MPI_Test(&req, &flag, MPI_STATUS_IGNORE);
        if (flag) {
            i->second.reset();
            m_mpiObjectMessages.erase(i++);
        } else {
            ++i;
        }
    }
    for (auto i = m_mpiMessages.begin(); i != m_mpiMessages.end();) {
        MPI_Request req = i->first;
        int flag;
        MPI_Test(&req, &flag, MPI_STATUS_IGNORE);
        if (flag) {
            delete i->second;
            m_mpiMessages.erase(i++);
        } else {
            ++i;
        }
    }
    if (m_shutdown) {
        return false;
    }
    return true;
}

bool Manager::checkMessages() {
    if (m_shutdown)
        return false;
    checkSends();
    int flag = 1;
    while (flag) {
        MPI_Status status;
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, m_comm, &flag, &status);
        if (flag) {
            switch (status.MPI_TAG) {
                case MPIRPC_TAG_SHUTDOWN:
                    m_shutdown = true;
                    handleShutdown();
                    checkSends();
                    return false;
                case MPIRPC_TAG_NEW:
                    registerRemoteObject();
                    break;
                case MPIRPC_TAG_INVOKE:
                    receivedInvocationCommand(std::move(status));
                    break;
                case MPIRPC_TAG_INVOKE_MEMBER:
                    receivedMemberInvocationCommand(std::move(status));
                    break;
                case MPIRPC_TAG_RETURN:
                    return true;
                default:
                    UserMessageHandler func = m_userMessageHandlers.at(status.MPI_TAG);
                    func(std::move(status));
            }
        }
    }
    return true;
}

void Manager::registerRemoteObject()
{
    ObjectInfo info;
    MPI_Status status;
    MPI_Recv(&info, 1, MpiObjectInfo, MPI_ANY_SOURCE, MPIRPC_TAG_NEW, m_comm, &status);
    registerRemoteObject(status.MPI_SOURCE, info.type, info.id);
}

void Manager::sync() {
    while (queueSize() > 0) { checkMessages(); } //block until this rank's queue is processed
    MPI_Request req;
    int flag;
    MPI_Ibarrier(m_comm, &req);
    do
    {
        MPI_Test(&req, &flag, MPI_STATUS_IGNORE);
        checkMessages();
    } while (!flag); //wait until all other ranks queues have been processed
}

void Manager::registerRemoteObject(int rank, TypeId type, ObjectId id) {
    ObjectWrapper<void> *a = new ObjectWrapper<void>();
    a->m_id = id;
    a->m_type = type;
    a->m_rank = rank;
    m_registeredObjects.push_back(a);
}

void Manager::shutdown() {
    int buf = 0;
    for (int i = 0; i < m_numProcs; ++i)
    {
        if (i != m_rank) {
            MPI_Bsend((void*) &buf, 1, MPI_INT, i, MPIRPC_TAG_SHUTDOWN, m_comm);
        }
    }
    m_shutdown = true;
    checkSends();
}

void Manager::handleShutdown()
{
    int buf;
    MPI_Status status;
    MPI_Recv(&buf, 1, MPI_INT, MPI_ANY_SOURCE, MPIRPC_TAG_SHUTDOWN, m_comm, &status);
    m_shutdown = true;
}

void Manager::receivedInvocationCommand(MPI_Status&& status)
{
    m_count++;
    int len;
    MPI_Get_count(&status, MPI_CHAR, &len);
    if (len != MPI_UNDEFINED) {
        std::vector<char>* buffer = new std::vector<char>(len);
        ParameterStream stream(buffer);
        MPI_Status recvStatus;
        MPI_Recv(stream.data(), len, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, m_comm, &recvStatus);
        FunctionHandle functionHandle;
        bool getReturn;
        stream >> functionHandle >> getReturn;
        FunctionBase *f = m_registeredFunctions[functionHandle];
        f->execute(stream, recvStatus.MPI_SOURCE, this, getReturn);
        delete buffer;
    }
}

void Manager::receivedMemberInvocationCommand(MPI_Status&& status) {
    m_count++;
    int len;
    MPI_Get_count(&status, MPI_CHAR, &len);
    if (len != MPI_UNDEFINED) {
        std::vector<char>* buffer = new std::vector<char>(len);
        ParameterStream stream(buffer);
        MPI_Status recvStatus;
        MPI_Recv(stream.data(), len, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, m_comm, &recvStatus);
        FunctionHandle functionHandle;
        ObjectId objectId;
        TypeId typeId;
        bool getReturn;
        stream >> typeId >> objectId >> functionHandle >> getReturn;
        FunctionBase *f = m_registeredFunctions[functionHandle];
        f->execute(stream, recvStatus.MPI_SOURCE, this, getReturn, getObjectWrapper(m_rank, typeId, objectId)->object());
        delete buffer;
    }
}

MPI_Comm Manager::comm() const
{
    return m_comm;
}

ObjectWrapperBase* Manager::getObjectOfType(mpirpc::TypeId typeId) const
{
    for (ObjectWrapperBase* i : m_registeredObjects)
    {
        if (i->type() == typeId)
            return i;
    }
    throw std::out_of_range("Object not found");
}

std::unordered_set<ObjectWrapperBase*> Manager::getObjectsOfType(TypeId typeId) const
{
    std::unordered_set<ObjectWrapperBase*> ret;
    for (ObjectWrapperBase* i : m_registeredObjects)
    {
        if (i->type() == typeId)
            ret.insert(i);
    }
    return ret;
}

ObjectWrapperBase* Manager::getObjectOfType(TypeId typeId, int rank) const
{
    for (ObjectWrapperBase* i : m_registeredObjects)
    {
        if (i->type() == typeId && i->rank() == rank)
            return i;
    }
    throw std::out_of_range("Object not found");
}

std::unordered_set< ObjectWrapperBase* > Manager::getObjectsOfType(TypeId typeId, int rank) const
{
    std::unordered_set<ObjectWrapperBase*> ret;
    for (ObjectWrapperBase* i : m_registeredObjects)
    {
        if (i->type() == typeId && i->rank() == rank)
            ret.insert(i);
    }
    return ret;
}

unsigned long long Manager::stats() const
{
    return m_count;
}

int Manager::numProcs() const
{
    return m_numProcs;
}

size_t Manager::queueSize() const
{
    return m_mpiObjectMessages.size() + m_mpiMessages.size();
}

ObjectWrapperBase* Manager::getObjectWrapper(int rank, TypeId tid, ObjectId oid) const {
    for (ObjectWrapperBase* i : m_registeredObjects)
        if (i->type() == tid && i->id() == oid && i->rank() == rank)
            return i;
    throw UnregisteredObjectException();
}

FunctionHandle Manager::FunctionBase::_idCounter = 0;

}
