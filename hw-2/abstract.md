HMWK 2 — From Wallet Seed to Ethereum Address: Hashing, Key Derivation and Cryptographic Primitives
Statistical properties

TASK
Check the practical cryptographic pipeline used in Ethereum wallets, from a mnemonic phrase to an Ethereum address and
briefly describe the role and purpose of the following standards and cryptographic primitives:

BIP-39 — mnemonic sentence and seed generation;
BIP-32 — hierarchical deterministic (HD) key derivation;
BIP-44 — standard derivation paths for hierarchical deterministic wallets;
secp256k1 — elliptic-curve public-key cryptography;
Keccak-256 — hash function used in Ethereum address derivation;
EIP-55 — mixed-case checksum encoding for Ethereum addresses.

In case try experiencing with software such as metamask, rabby, exodus, rainbow, ... or similar
to understand the role of each iece of technology ...

A simplified wallet pipeline can be represented as:

Mnemonic ↓ BIP-39 ↓ 
Seed  BIP-32 / BIP-44 ↓ 
Private Key ↓ secp256k1 ↓ 
Public Key ↓ Keccak-256 ↓ 
Ethereum Address ↓ EIP-55 ↓ 
Checksummed Address

For the purpose of this homework, distinguish clearly between:

key generation and key derivation;
public-key cryptography;
hashing;
checksums;
encoding.
Also briefly discuss the role of each primitive and how the different steps fit together
and the statistical properties.


Hints: TEMPORAL AND HISTORICAL CONTEXT

Briefly place these technologies in chronological context.

Discuss the approximate historical order of:
RSA — public-key cryptography, introduced in 1977;
SHA-2 / SHA-512 — standardized in 2001;
BIP-32 — hierarchical deterministic wallets, 2012;
BIP-39 — mnemonic code for generating deterministic wallets, 2013;
BIP-44 — multi-account hierarchical deterministic wallets, 2014;
SHA-3 — standardized by NIST in 2015, based on Keccak;
EIP-55 — Ethereum mixed-case checksum addresses, 2016.
Also explain the distinction between:

SHA-2 / SHA-512 and Keccak / SHA-3

They are related historically and both belong to the family of cryptographic hash functions, but Keccak is the underlying design selected for SHA-3, while Ethereum uses Keccak-256, which is not identical to the standardized SHA3-256 function (because of the padding/domain-separation convention.)

Why does the Ethereum wallet pipeline use several different cryptographic techniques instead of using one single algorithm for everything?
Explain briefly why:
BIP-39/BIP-32/BIP-44 → key management and derivation
secp256k1 → public/private key relationship and digital signatures
Keccak-256 → hashing
EIP-55 → checksum/error detection
serve different purposes in the overall system