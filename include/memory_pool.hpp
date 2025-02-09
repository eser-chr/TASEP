#include <stdexcept>

#pragma once

class MemoryPool {
public:
    MemoryPool(void* memory, std::size_t size)
        : memory_start_(static_cast<char*>(memory)), pool_size_(size), offset_(0)
    {
        if (!memory_start_ || pool_size_ == 0)
            throw std::invalid_argument("Invalid memory pool initialization");
    }

    // Allocate a chunk of memory from the pool.
    // Returns a pointer to the allocated block, or throws std::bad_alloc.
    void* allocate(std::size_t bytes) {
        if (offset_ + bytes > pool_size_)
            throw std::bad_alloc();
        void* ptr = memory_start_ + offset_;
        offset_ += bytes;
        return ptr;
    }

    // Optionally, add a reset() method if you want to reclaim the whole pool at once.
    void reset() { offset_ = 0; }

    // For debugging
    std::size_t used() const { return offset_; }
    std::size_t total_size() const { return pool_size_; }

private:
    char* memory_start_;
    std::size_t pool_size_;
    std::size_t offset_;
};


template <typename T>
class PoolAllocator {
public:
    using value_type = T;

    // Constructor: store a pointer to the shared MemoryPool.
    PoolAllocator(MemoryPool* pool = nullptr) : pool_(pool) {}

    // Copy constructor template for rebinding.
    template <typename U>
    PoolAllocator(const PoolAllocator<U>& other) : pool_(other.pool()) {}

    T* allocate(std::size_t n) {
        if (!pool_)
            throw std::bad_alloc();
        // Calculate bytes needed and use the pool to allocate.
        std::size_t bytes = n * sizeof(T);
        return static_cast<T*>(pool_->allocate(bytes));
    }

    void deallocate(T* p, std::size_t n) {
        // For a bump allocator, individual deallocation is typically a no-op.
        // More complex pool management might track and reuse freed memory.
    }

    MemoryPool* pool() const { return pool_; }

    // Allocators are considered equal if they refer to the same pool.
    template <typename U>
    bool operator==(const PoolAllocator<U>& other) const {
        return pool_ == other.pool();
    }
    template <typename U>
    bool operator!=(const PoolAllocator<U>& other) const {
        return !(*this == other);
    }

private:
    MemoryPool* pool_;

    // Allow access to pool_ for allocators of other types.
    template <typename U>
    friend class PoolAllocator;
};

