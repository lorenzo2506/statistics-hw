HMWK 1 — Example App: Public-Key Cryptography PRIVATE KEY → PUBLIC KEY in Ethereum

Notes: Elliptic-curve scalar multiplication can be viewed as a cryptographic one-way function.
It is deterministic and easy to compute in the forward direction: private key k → public key Q = kG
while recovering k from Q is computationally difficult. It therefore shares some properties with cryptographic 
hash functions, but it is based on the Elliptic Curve Discrete Logarithm Problem (ECDLP) and is used for public-key 
cryptography and digital signatures.

The goal is to understand the basic idea of elliptic-curve scalar multiplication, not to implement a real cryptographic system.

TASK
Write a small JavaScript/html/css program showing how a private key can be used to generate a public key using elliptic-curve 
scalar multiplication. (In Ethereum, the relevant elliptic curve is secp256k1).

Hints: For this exercise, use a small toy curve: y² = x³ + 7 (mod 17)
For example, use for instance the following base point: G = (1, 5) Choose for instance a small private key: k = 6
The public key is: Q = kG  Define the toy curve and the base point G. Define the private key k.
Implement the basic point-addition operation. Compute and print: G, 2G, 3G, ..., kG. Or use  afaster method for a larger k.
Print the final point as the TOY PUBLIC KEY.

OPTIONAL-GURU LEVEL :-)  Make a slown-down visual representation of the algo (may use the canvas or SVG, whatever you like better)

== STATISTICAL EXPLORATION ==

Repeat the experiment with different values of k and, if useful, different valid base points.
Explore simple statistical properties of the resulting points, for example:

distribution of x and y coordinates;
frequency of individual values;
distribution of bits;
apparent uniformity;
correlations between input and output;
possible avalanche-like behaviour.
Discuss which of these properties are meaningful for a cryptographic one-way function and which are properties specifically associated with cryptographic hash functions.

---------------------------------------------------------
