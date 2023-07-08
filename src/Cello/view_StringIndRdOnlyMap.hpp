// See LICENSE_CELLO file for license and copyright information

/// @file     view_StringIndRdOnlyMap.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 30 2019
/// @brief    Declaration and implementation of the StringIndRdOnlyMap class

#include <array>
#include <functional> // hash
#include <limits>
#include <memory>
#include <string>
#include <vector>

#ifndef VIEW_STRING_IND_RD_ONLY_MAP_HPP
#define VIEW_STRING_IND_RD_ONLY_MAP_HPP

//----------------------------------------------------------------------

/// @def    INVERSE_LOAD_FACTOR
/// @brief  the load factor specifies the fraction of the capcity of the Hash
///         table that is filled. This should be an integer.
///
/// This should not be 1. Generally, the larger this is, the fewer collisions
/// there are, but the more memory is required. If this is too large, it may
/// actually slow down lookups. 
#define INVERSE_LOAD_FACTOR 2

//----------------------------------------------------------------------

class StringIndRdOnlyMap{

  /// @class    StringIndRdOnlyMap
  /// @ingroup  View
  /// @brief    [\ref View] The primary purpose of this is to associate string
  ///           keys with indices. If you have `n` keys, then they will be
  ///           associated with an ordered set values from `0` to `n-1`.
  ///
  /// This is primarily intended to be used in the implementation of the array
  /// map (but may be broadly useful for other applications).
  ///
  /// To make instances cheap to copy, we store the actual hash table in a
  /// shared_ptr. There is no downside to doing this since the hash table
  /// can never be mutated after it is constructed.
  ///
  /// We implement our own simple hash-table instead of using one of the standard
  /// library datastructures (i.e. ``std::map``) to avoid an extra level of
  /// indirection. This uses a simple open-addressed hash-table (with linear
  /// probing to address collisions). Since we don't allow the contents to be
  /// mutated, most of the negatives of an open-addressed hash-table do not
  /// apply.
  ///
  /// We prioritized making the implementation fast for looking up the index
  /// associated with a given key. This also supports looking up the value of
  /// the `i`th key, with the `key` method, (basically a reverse lookup). But,
  /// the `key` method is mostly for convenience and is slow. If we ever
  /// the `key` method to be faster, this class could additionally track a
  /// shared read-only array that stores the order of keys.
  ///
  /// There is definitely room for optimizing this implementation:
  /// - We could be smarter about the order that we insert keys into the table
  ///   (in the constructor) to minimize the search time.
  /// - We might be able to come up with a better hash function
  /// - If we limit the max key size we could have `KVpair` store keys in-place.
  ///   For example if the max key size never exceeds ~22 characters and
  ///   there are never more than ~128 entries, then we could define `KVpair`
  ///   as `KVpair{ char key[23]; uint8_t val; };`. (This reduces indirection).
  ///   If `sizeof(KVpair)` is `24` or `32`, we also could fit 3 or 2 entries
  ///   per cache line.

public:
  typedef uint16_t val_t;
  struct KVpair{ std::string key; val_t val; };

public:
  /// Constructor
  StringIndRdOnlyMap()
    : data_(nullptr),
      num_keys_(0),
      capacity_(0),
      hash_()
  {}

  /// Constructor
  inline StringIndRdOnlyMap(const std::vector<std::string> &keys) noexcept;

  /// returns the value associated with the key
  ///
  /// The program aborts if the key isn't in the table.
  inline val_t operator[](const std::string& key) const noexcept
  { 
    const KVpair& tmp = find_match_or_empty_(key);
    if (tmp.key == ""){
      ERROR1("StringIndRdOnlyMap::operator[]",
             "there is no key called %s", key.c_str());
    }
    return tmp.val;
  }

  /// returns the value associated with the key
  ///
  /// The program aborts if the key isn't in the table.
  inline val_t at(const std::string& key) const noexcept
  { return (*this)[key]; }

  /// returns the number of entries in the map
  inline std::size_t size() const noexcept {return num_keys_;}

  /// checks if the map contains a key
  inline bool contains(const std::string& key) const noexcept
  { return (find_match_or_empty_(key).key == key); }

  /// Return the ith key (this is effectively a reverse lookup)
  inline const std::string& key(val_t i) const noexcept{
    if (i >= (num_keys_)){
      ERROR2("StringIndRdOnlyMap::key",
             "Can't find key number %d. There only %d keys.",
             (int) i, (int)num_keys_); 
    }

    /// iterate over over the entire capacity of the hash table
    for (uint32_t j = 0; j < capacity_; j++){
      const KVpair& pair = data_.get()[j];
      if ((pair.key != "") & (pair.val == i)){ return pair.key; }
    }

    // the following should be unreachable
    ERROR2("StringIndRdOnlyMap::find_match_or_empty_",
           "something went horribly wrong while searching for key number %d. "
           "the max number of keys is %d",
           (int) i, (int)num_keys_);
  }

private:

  const KVpair& find_match_or_empty_(const std::string& key) const {
    std::size_t hash_code = hash_(key) % capacity_;

    for (std::size_t i = hash_code; i < capacity_; i++){
      const KVpair& pair = data_.get()[i];
      if ((pair.key == key) | (pair.key == "")) {return pair;}
    }

    for (std::size_t i = 0; i < hash_code; i++){
      const KVpair& pair = data_.get()[i];
      if ((pair.key == key) | (pair.key == "")) {return pair;}
    }

    ERROR("StringIndRdOnlyMap::find_match_or_empty_",
          "something went horribly wrong - this should be unreachable");
  }

private:
  std::shared_ptr<KVpair> data_;
  std::uint32_t num_keys_;
  std::uint32_t capacity_;
  std::hash<std::string> hash_;
};

//----------------------------------------------------------------------

inline StringIndRdOnlyMap::StringIndRdOnlyMap
(const std::vector<std::string> &keys) noexcept
{

  if ((std::size_t)(std::numeric_limits<val_t>::max()) < keys.size()){
    ERROR("StringIndRdOnlyMap::StringIndRdOnlyMap",
          "There are WAY too many keys");
  }

  // determine capacity of hash table
  capacity_ = 0;
  uint32_t target_capacity = (uint32_t) (keys.size() * INVERSE_LOAD_FACTOR);

  // we want the capacity to be a prime number
  std::array<uint32_t,12> allowed_caps = {7,  19,  31,  41,
                                          53, 61,  71,  83,
                                          89, 101, 113, 127};

  if (target_capacity == 0){
    ERROR("StringIndRdOnlyMap::StringIndRdOnlyMap",
          "can't construct a map without any keys");
  }
  
  for (const uint32_t &val : allowed_caps){
    if (target_capacity < val){
      capacity_ = val;
      break;
    }
  }
  if (capacity_ == 0){
    ERROR("StringIndRdOnlyMap::StringIndRdOnlyMap",
          "There are too many keys for the current implementation");
  }
  num_keys_ = keys.size();

  // allocate hash table (all KVpairs should be default constructed)
  KVpair* ptr = new KVpair[capacity_]();
  if (ptr == nullptr){
    ERROR("StringIndRdOnlyMap::StringIndRdOnlyMap",
          "Something went wrong while allocating data");
  }

  data_ = std::shared_ptr<KVpair>(ptr, std::default_delete<KVpair[]>());

  // fill in the hash table:
  for (std::size_t i = 0; i < keys.size(); i++){
    const std::string& key = keys[i];
    if (key == ""){
      ERROR("StringIndRdOnlyMap::StringIndRdOnlyMap",
            "Empty strings are not allowed to be keys");
    }
    KVpair& pair = const_cast<KVpair&>(find_match_or_empty_(key));
    if (pair.key == key){
      ERROR1("StringIndRdOnlyMap::StringIndRdOnlyMap",
             "There can't be more than 1 key called %s.",
             key.c_str());
    }
    pair = {key, (val_t)i};
  }
}


#endif /* VIEW_STRING_IND_RD_ONLY_MAP_HPP */
