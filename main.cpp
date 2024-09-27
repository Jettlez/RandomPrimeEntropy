#include <bitset>
#include <iostream>
#include <vector>
#include <windows.h>
#include <thread>
#include <chrono>
#include <unordered_set>


/**
 * The size of each chunk for bitset, defined as 1 MB (1024 * 1024 bytes).
 *
 * This constant is used to partition the bitset storage for managing a large
 * amount of prime number calculations efficiently.
 */
constexpr size_t BITSET_CHUNK_SIZE = 1024 * 1024;

/**
 * Returns the width of the primary display monitor in pixels.
 *
 * @return int Width of the primary display monitor.
 */
constexpr int GetScreenWidth() {
    return GetSystemMetrics(SM_CXSCREEN);
}

/**
 * @brief Retrieves the height of the primary display screen.
 *
 * This function returns the height of the screen in pixels by
 * calling the GetSystemMetrics function with the SM_CYSCREEN flag.
 * The value is determined at compile time since the function is
 * marked as constexpr.
 *
 * @return The height of the primary display screen in pixels.
 */
constexpr int GetScreenHeight() {
    return GetSystemMetrics(SM_CYSCREEN);
}

/**
 * @brief Retrieves the maximum number of hardware threads supported by the system.
 *
 * This function uses the `std::thread::hardware_concurrency` method to determine
 * the maximum number of concurrent threads supported by the system hardware.
 * The value returned is an approximation and may not necessarily equal the number
 * of threads the system can efficiently handle.
 *
 * @return The maximum number of concurrent threads supported by the hardware.
 */
constexpr unsigned int GetMaxThreads() {
    return std::thread::hardware_concurrency();
}

/**
 * Generates an unsigned long integer based on screen resolutions, current time, and system tick count.
 *
 * This function computes a unique value using the screen width, current time, screen height,
 * and system tick count. It adds together these values to generate a somewhat random and unique result.
 *
 * @return An unsigned long integer combining screen resolutions, current time, and system tick count.
 */
constexpr unsigned long long EntropyGen() {
    return (GetScreenWidth() + time(nullptr) + GetScreenHeight() + time(nullptr) + GetTickCount64());
    //returns the constexpr resolutions multiplied by each other plus the TickCount of the computer
}

/**
 * @brief Marks the non-prime numbers in a given segment as false in the passed boolean vector.
 *
 * This function implements the segmented Sieve of Eratosthenes algorithm. For a given range
 * [start, end], it marks all non-prime numbers within this range as false in the provided
 * is_prime vector. The vector is assumed to have already been initialized such that its
 * size is at least `end + 1` and all values are initially set to true.
 *
 * @param start The max's starting number of the range to be processed.
 * @param end The max's ending number of the range to be processed.
 * @param is_prime Reference to a vector of booleans representing the primality of numbers.
 *                 If is_prime[i] is true, then 'i' is assumed to be a prime number.
 */
void sieve_segment(unsigned int start, unsigned int end, std::vector<bool> &is_prime) {
    //takes a starting number , ending number and a reference to a vector that contains a list of whether numbers are prime or not
    for (unsigned int num = start; num <= end; num++) {
        //loops every number in the specified range from start to end
        if (is_prime[num]) {
            // checks if the number is marked with 'a' is_prime value
            for (unsigned int multiple = num * num; multiple < is_prime.size(); multiple += num) {
                //starts at num*num since any smaller multiple of num would be marked prime, multiple is increased by num each iteration marking each multiple of num
                is_prime[multiple] = false;
                //for each multiple if is_prime[multiple] is false that number is not a prime
            }
        }
    }
}

/**
 * @brief Counts the number of prime and non-prime numbers in a specified range.
 *
 * This function takes a start value and an end value, along with a vector indicating
 * if a number is prime, then counts how many numbers in the range [start, end) are prime
 * and how many are non-prime. The results are stored in `local_count` for prime
 * numbers and `non_prime_count` for non-prime numbers.
 *
 * @param start The start value of the range.
 * @param end The end value of the range (exclusive).
 * @param is_prime A vector<bool> that indicates if the number is prime.
 * @param local_count Reference to an unsigned int that will hold the count of prime numbers.
 * @param non_prime_count Reference to an unsigned int that will hold the count of non-prime numbers.
 */
void count_primes_and_non_primes(unsigned int start, unsigned int end, const std::vector<bool> &is_prime,
                                 unsigned int &local_count, unsigned int &non_prime_count) {
    //takes a start value and an end value, and a reference to the is_prime vector, as well as a reference to an unsigned int that keeps count of prime numbers in the range given
    for (unsigned int i = start; i < end; i++) {
        //loops over each number from start to end-1
        if (is_prime[i]) {
            //if the number it check is_prime then the number is prime
            local_count++; //+1 to the local_count to track number of prime numbers found in the range
        } else {
            // else if prime is not true then non_prime is discovered and 1 is added to the count
            non_prime_count++;
        }
    }
}

/**
 * @class RandomNumberGenerator
 * @brief A simple random number generator using a linear congruential generator (LCG) algorithm.
 *
 * This class utilizes the linear congruential generator algorithm to produce a sequence of
 * pseudo-random numbers.
 */
class RandomNumberGenerator {
    /**
     * @brief Seed value for the random number generator.
     *
     * Used as the starting point for generating pseudo-random numbers.
     */
private:
    unsigned long long seed;
    const unsigned long long m = 4294967296;
    /**
     * @brief Multiplier parameter for the linear congruential generator (LCG).
     *
     * This constant represents the multiplier used in the LCG algorithm for
     * generating pseudo-random numbers. The value is chosen to meet certain
     * criteria for randomness.
     */
    const unsigned long long a = 6633142268346709350ULL;
    /**
     * @brief Increment value used in the random number generation algorithm.
     *
     * This value is used to ensure the pseudo-random sequence generated by the
     * RandomNumberGenerator class is sufficiently varied. It is part of the
     * linear congruential generator (LCG) algorithm parameters, specifically
     * representing the increment value (c) in the equation:
     *
     *     X_(n+1) = (a * X_n + c) % m
     */
    const unsigned long long c = 1ULL << 63; // Increment

    /**
     * @brief Initializes the RandomNumberGenerator with the given seed value.
     *
     * This constructor sets the internal seed for the linear congruential generator (LCG)
     * algorithm used in generating random numbers.
     *
     * @param seed The initial seed value used in generating random numbers.
     */
public:
    explicit RandomNumberGenerator(unsigned long long seed) : seed(seed) {
    }

    /**
     * @brief Generates a random number based on the current seed.
     *
     * This function uses a linear congruential generator (LCG) to produce a pseudo-random number.
     * The formula used is: seed = (a * seed + c) % m, where:
     * - 'a' is a constant multiplier.
     * - 'c' is an increment.
     * - 'm' is the modulus (2^32 in this case).
     *
     * @return The next pseudo-random number in the sequence.
     */
    unsigned long long getRandomNumber() {
        seed = (a * seed + c) % m;
        return seed;
    }
};

/**
 * Tests the provided random number generator for duplicate detection and performance.
 *
 * This function generates a specified number of random numbers using the provided
 * random number generator and checks for any duplicates.
 *
 * @param rng Reference to an instance of a RandomNumberGenerator.
 * @param num_samples Number of random numbers to generate and test for duplicates.
 */
void testRandomNumberGenerator(RandomNumberGenerator &rng, int num_samples) {
    std::unordered_set<unsigned long long> seen;

    // Generate numbers and check for duplicates
    for (int i = 0; i < num_samples; ++i) {//iterate from 0 to num_samples generating a new random number each time
        unsigned long long num = rng.getRandomNumber(); //random number generated by calling getRandomNumber from the rng instance then place it into num
        if (seen.contains(num)) { // "seen" our unordered_set has previously generated numbers; contain checks if the number is already in "seen"
            std::cout << "Duplicate found: " << num << std::endl; //if a duplicate is found the program addresses it, which number then breaks the program
            break; // Stop on first duplicate
        }
        seen.insert(num); //if the number wasn't found in "seen" it is added to "seen" to track in the future.
    }
}

/**
 * @class PrimeNumberGenerator
 * @brief Class for generating prime numbers up to a specified maximum using the Sieve of Eratosthenes algorithm.
 *
 * This class generates prime numbers up to a given maximum value and calculates the total number of prime
 * and non-prime numbers within that range. The process is optimized by using multiple threads to count
 * prime numbers in parallel.
 *
 * @param max The maximum value up to which prime numbers are to be generated.
 */
class PrimeNumberGenerator {
    /**
     * Constructor that initializes the PrimeNumberGenerator object.
     *
     * This function performs the following operations:
     * - Starts a high-resolution timer to measure the execution time.
     * - Initializes a bitset vector to mark prime numbers up to the given maximum value
     *   using the Sieve of Eratosthenes algorithm.
     * - Implements the Sieve of Eratosthenes to mark non-prime numbers.
     * - Divides the range of numbers into chunks and counts prime and non-prime numbers
     *   using multithreading for parallel processing.
     * - Waits for all threads to finish execution and then sums the local counts of prime
     *   and non-prime numbers.
     * - Stops the timer and calculates the duration taken for the entire computation.
     * - Displays the total number of primes, total number of non-primes,
     *   and the execution time.
     * - Initializes a RandomNumberGenerator object with the total number of primes and
     *   an entropy generator.
     * - Tests the RandomNumberGenerator with a specified number of random numbers.
     *
     * @param max The maximum value up to which prime numbers are to be generated.
     */
public:
    explicit PrimeNumberGenerator(unsigned int max) {
        auto start_time = std::chrono::high_resolution_clock::now(); // Start timer

        size_t total_chunks = (max + BITSET_CHUNK_SIZE - 1) / BITSET_CHUNK_SIZE;
        std::vector<std::bitset<BITSET_CHUNK_SIZE> > is_prime(total_chunks);

        // Initialize is_prime
        for (size_t i = 0; i < total_chunks; ++i) {
            is_prime[i].set(); // Set all bits to true
        }
        is_prime[0].reset(0); // 0 is not prime
        is_prime[0].reset(1); // 1 is not prime

        // Implement the Sieve of Eratosthenes for the first sqrt(max)
        for (unsigned int i = 2; i * i <= max; i++) {
            //loops over potential prime numbers starting from 2, i * i <= max mean it only checks for factors up to the square root of max
            if (is_prime[i / BITSET_CHUNK_SIZE].test(i % BITSET_CHUNK_SIZE)) {
                //if the number it checked is marked in the prime vector than 'i' is prime
                for (unsigned int j = i * i; j <= max; j += i) {
                    //the loops start from i * i which is the first multiple of i that is non-prime, smaller multiples like 2i would already be marked by smaller prime factors,
                    is_prime[j / BITSET_CHUNK_SIZE].reset(j % BITSET_CHUNK_SIZE);
                    //multiples of i are marked as non-prime since all prime values multiplied by themselves are non-prime (5*5 = 10)
                }
            }
        }

        // Count total primes
        unsigned int total_primes = 0; //total number of prime numbers found
        unsigned int total_non_primes = 0; //total number of non_primes found
        const unsigned int num_threads = GetMaxThreads();
        //number of concurrent threads the cpu can support
        std::vector<std::thread> threads; //store the thread objects
        std::vector<unsigned int> local_counts(num_threads, 0); //keeps track of prime counts for each thread
        std::vector<unsigned int> non_prime_counts(num_threads, 0); //keeps track of non_prime counts for each thread

        unsigned int range_per_thread = max / num_threads;

        for (unsigned int i = 0; i < num_threads; i++) {
            //loops 0 to num_threads - 1, creating a seperated thread.
            unsigned int start = i * range_per_thread + 2;
            //calculates the starting index for the thread; then each thread is given a range based on that segment
            unsigned int end = (i == num_threads - 1) ? max + 1 : start + range_per_thread;
            // determines where the thread stops; with the last thread being max + 1 to make sure it gets all the numbers

            threads.emplace_back([start, end, &is_prime, &local_counts, &non_prime_counts, i]() {
                //emplace_back is used to create a new thread by passing a lambda function; the lambda function gets a start, end, and the is_prime bitset as well as the local_counts and non_prime_counts vectors, then the threads index 'i'
                unsigned int local_count = 0; //thread is given its own individual local_count
                unsigned int non_prime_count = 0; //thread is given its own individual non_prime_count
                for (unsigned int j = start; j < end; j++) {
                    //loop through the range given to the thread
                    if (is_prime[j / BITSET_CHUNK_SIZE].test(j % BITSET_CHUNK_SIZE)) {
                        //access the bitset for the number j; assumes you have BITSET_CHUNK_SIZE space.
                        local_count++;
                    } else {
                        non_prime_count++;
                    }
                }
                local_counts[i] = local_count; //each thread stores its local_counts vector into its thread ID
                non_prime_counts[i] = non_prime_count;
                //each thread stores its non_prime_counts vector into its thread ID
            });
        }

        for (auto &thread: threads) {
            //we wait for each thread to finish its job (specifically the last thread) by using join()
            thread.join();
        }

        for (const auto &count: local_counts) {
            total_primes += count; //after the threads finish combine all the total_primes count
        }

        for (const auto &count: non_prime_counts) {
            total_non_primes += count; //after the threads finish combine all the total_primes count
        }

        auto end_time = std::chrono::high_resolution_clock::now(); // End timer
        std::chrono::duration<double> duration = end_time - start_time; // Calculate speed

        std::cout << "Total primes found: " << total_primes << std::endl;
        std::cout << "Total non primes found: " << total_non_primes << std::endl;
        std::cout << "Time taken: " << duration.count() << " seconds" << std::endl;
        // Display how long the program took to count all threads numbers.

        // Initialize RandomNumberGenerator with total_primes and increment from EntropyGen
        RandomNumberGenerator rng(EntropyGen());
        //initialize RandomNumberGenerator and supply it the total primes, the Entropy Gen, as well as a mix of both

        // Generate random numbers
        testRandomNumberGenerator(rng, 1000000000);
    }
};

/**
 * @brief The main entry point for the program.
 *
 * This function generates entropy and uses it to initialize a PrimeNumberGenerator
 * with a maximum value. The entropy is divided by 16.
 * The PrimeNumberGenerator is then used to calculate prime numbers up to the maximum value.
 *
 * @return int Returns 0 upon successful execution.
 */
int main() {
    PrimeNumberGenerator((EntropyGen() / 16));
    // set the max amount of number you want, or set it to EntropyGen or really any function that returns a large int
    std::cout << "Press Enter to exit...";
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear any leftover input
    std::cin.get(); // Wait for the user to press Enter
    return 0;
}
