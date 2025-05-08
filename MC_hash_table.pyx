#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 23 20:54:39 2025

@author: mattc
"""
from libc.string cimport memcmp
from libc.stdlib cimport malloc, calloc, free
from libc.stdint cimport uint64_t, uint8_t


cdef struct Kmer128:
    uint64_t high
    uint64_t low

cdef inline uint8_t base_to_bits(char base):
    if base == b'A':
        return 0
    elif base == b'C':
        return 1
    elif base == b'G':
        return 2
    elif base == b'T':
        return 3
    else:
        return 0  # or raise error

cdef Kmer128 encode_69mer(const char* s):
    cdef Kmer128 k
    cdef int i
    k.low = 0
    k.high = 0
    for i in range(0, 32):
        k.low = (k.low << 2) | base_to_bits(s[i])
    for i in range(32, 64):
        k.high = (k.high << 2) | base_to_bits(s[i])
    for i in range(64, 69):
        k.high = (k.high << 2) | base_to_bits(s[i])
    return k

cdef struct Entry:
    Kmer128 key
    int value
    Entry* next
    
cdef class DNAHashTable:
    cdef Entry** table
    cdef int size

    def __cinit__(self, int size=1024):
        self.size = size
        self.table = <Entry**>calloc(size, sizeof(Entry*))

    def __dealloc__(self):
        cdef int i
        cdef Entry* curr
        cdef Entry* tmp
        for i in range(self.size):
            curr = self.table[i]
            while curr:
                tmp = curr.next
                free(curr)
                curr = tmp
        free(self.table)

    cpdef void insert_69mer(self, str kmer, int value):
        cdef bytes kmer_bytes = kmer.encode("ascii")
        cdef const char* c_kmer = kmer_bytes
        cdef Kmer128 kmer_key = encode_69mer(c_kmer)
        cdef unsigned int idx = (kmer_key.low ^ kmer_key.high) % self.size

        cdef Entry* new_entry = <Entry*>malloc(sizeof(Entry))
        new_entry.key = kmer_key
        new_entry.value = value
        new_entry.next = self.table[idx]
        self.table[idx] = new_entry

    cpdef int lookup_69mer(self, str kmer):
        cdef bytes kmer_bytes = kmer.encode("ascii")
        cdef const char* c_kmer = kmer_bytes
        cdef Kmer128 kmer_key = encode_69mer(c_kmer)
        cdef unsigned int idx = (kmer_key.low ^ kmer_key.high) % self.size

        cdef Entry* curr = self.table[idx]
        while curr:
            if curr.key.low == kmer_key.low and curr.key.high == kmer_key.high:
                return curr.value
            curr = curr.next
        return -1
