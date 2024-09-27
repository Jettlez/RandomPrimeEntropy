# RandomPrimeEntropy
<u><i>"Put your hand on a hot stove for a minute"</i></u> - Evil Einstein
# How it works?
- first entropy is gathered from the system using windows functions like getsystemmetrics() and gettickcount64() into a constexpr called EntropyGen()
- using the <u>Sieve of Eratosthenes</u> we can then loop over the numbers generated from EntropyGen and finds all prime and non prime numbers and places them into individual counts
- using the values gathered from the prime functions along with EntropyGen(), we can multiply, add, divide, square, really do whatever we want to generate seeds and use those seeds to generate values.
- Win
