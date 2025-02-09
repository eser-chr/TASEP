#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <new>
#include <stdexcept>
#include <vector>

// A simple bump-allocator that uses a fixed block of memory.
template <typename T>
class BumpAllocator {
   public:
    using value_type = T;

    // Pointer types (C++11 compliant)
    using pointer = T*;
    using const_pointer = const T*;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;

    // Rebind allocator to type U
    template <typename U>
    struct rebind {
        using other = BumpAllocator<U>;
    };

    // Constructor: pass a pointer to pre-allocated memory and its size (in bytes)
    BumpAllocator(void* memory, size_type memory_size)
        : memory_start_(static_cast<char*>(memory)), memory_size_(memory_size), offset_(0) {
        // Ensure that the memory block is not null if memory_size > 0.
        if (memory_size_ == 0) {
            throw std::invalid_argument("Memory size must be greater than zero.");
        }
    }

    // Copy constructor: share the same memory block for simplicity.
    template <typename U>
    BumpAllocator(const BumpAllocator<U>& other)
        : memory_start_(other.memory_start_),
          memory_size_(other.memory_size_),
          offset_(other.offset_) {}

    // Allocate memory for n objects of type T.
    pointer allocate(size_type n, const void* hint = 0) {
        size_type bytes_needed = n * sizeof(T);
        if (offset_ + bytes_needed > memory_size_) {
            throw std::bad_alloc();
        }
        pointer result = reinterpret_cast<pointer>(memory_start_ + offset_);
        offset_ += bytes_needed;
        // For simplicity, we don't call any constructors here.
        return result;
    }

    // Deallocate memory: this is a no-op in this bump allocator.
    void deallocate(pointer p, size_type n) {
        // In a bump allocator, individual deallocation is not supported.
        // For a real implementation, you might want to add logic to reuse freed memory.
    }

    // Equality operators (allocators with the same memory block are considered equal)
    template <typename U>
    bool operator==(const BumpAllocator<U>& other) const {
        return memory_start_ == other.memory_start_;
    }

    template <typename U>
    bool operator!=(const BumpAllocator<U>& other) const {
        return !(*this == other);
    }

    // Reset the allocator (e.g., when you're done with the container, you can reclaim all memory)
    void reset() { offset_ = 0; }

    // Accessors for debugging
    char* memory_start() const { return memory_start_; }
    size_type memory_size() const { return memory_size_; }
    size_type used() const { return offset_; }

   private:
    // Pointer to the beginning of the memory block.
    char* memory_start_;
    // Total size of the memory block in bytes.
    size_type memory_size_;
    // Current offset in the block.
    size_type offset_;

    // Grant access to other instantiations of the template.
    template <typename U>
    friend class BumpAllocator;
};

int main() {
    // Allocate a block of memory (for example, 1024 bytes)
    const std::size_t block_size = 1024;
    void* memory_block = std::malloc(block_size);
    if (!memory_block) {
        std::cerr << "Failed to allocate memory block." << std::endl;
        return 1;
    }

    try {
        // Create a vector of ints that uses the bump allocator.
        // Note: It's important to ensure the memory block is large enough for the vector's needs.
        std::vector<int, BumpAllocator<int>> vec(BumpAllocator<int>(memory_block, block_size));

        // Reserve some space if you know the approximate number of elements:
        vec.reserve(50);

        // Fill the vector
        for (int i = 0; i < 50; ++i) {
            vec.push_back(i);
        }

        // Use the vector
        for (const auto& elem : vec) {
            std::cout << elem << " ";
        }
        std::cout << "\nUsed memory: " << vec.get_allocator().used() << " out of "
                  << vec.get_allocator().memory_size() << " bytes.\n";
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << "\n";
        std::free(memory_block);
        return 1;
    }

    // Free the pre-allocated memory block
    std::free(memory_block);

    return 0;
}
