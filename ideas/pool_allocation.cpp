#include <cstddef>
#include <stdexcept>
#include <cstdlib>
#include <cstdint>
#include <memory>
#include <new>
#include <vector>
#include <iostream>

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



template <typename T>
class BasicIteration /*: public AbstractIteration<T>*/ {
public:
    BasicIteration(MemoryPool* pool)
      : DATA(PoolAllocator<uint8_t>(pool)),
        TIMES(PoolAllocator<T>(pool)),
        ACTION(PoolAllocator<uint8_t>(pool)),
        SIDE(PoolAllocator<uint16_t>(pool))
    {
        // Now all vectors use the same underlying memory pool.
    }

    // Example vectors with custom allocators.
    std::vector<uint8_t, PoolAllocator<uint8_t>> DATA;
    std::vector<T, PoolAllocator<T>> TIMES;
    std::vector<uint8_t, PoolAllocator<uint8_t>> ACTION;
    std::vector<uint16_t, PoolAllocator<uint16_t>> SIDE;

private:
    void append_trajectory() {
        // Implementation that uses the vectors.
    }
};

int main() {
    // Pre-allocate a block of memory (for example, 1024 bytes or more as needed)
    const std::size_t poolSize = 4096;  // adjust according to your needs
    void* memoryBlock = std::malloc(poolSize);
    if (!memoryBlock) {
        std::cerr << "Failed to allocate memory pool.\n";
        return 1;
    }

    // Create the memory pool.
    MemoryPool pool(memoryBlock, poolSize);

    // Now create an instance of BasicIteration with the shared pool.
    BasicIteration<double> iteration(&pool);

    // Use the vectors as usual.
    iteration.DATA.push_back(42);
    iteration.TIMES.push_back(3.14);
    iteration.ACTION.push_back(255);
    iteration.SIDE.push_back(1000);

    // Optionally, print debug info.
    std::cout << "Memory used: " << pool.used() << " out of " << pool.total_size() << "\n";

    std::free(memoryBlock);
    return 0;
}