// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoFieldAdaptor.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2022-04-01
/// @brief    [\ref Enzo] Declaration of the EnzoFieldAdaptor class
///
/// The `EnzoFieldAdaptor` class is used to provide a common interface for
/// accessing field-like arrays from the `Field` class and `EnzoEFltArrayMap`.
/// We might be able to do away with this class once we make it easier to
/// construct an `EnzoEFltArrayMap` instance from a `Field` instance.
///
/// The primary motivation for this class is flexibility. Previously, reusable
/// procedures for computing derived fluid quantities (e.g. EnzoComputePressure)
/// accepted a `Block` instance as the primary argument and operated upon the
/// associated data. While this approach works, it's somewhat limiting since
/// `Block` objects carry a lot of unnecessary information. There are
/// situations where we would like to be able to compute a quantity without
/// needing the data to be directly associated with a `Block`.
///
/// Consider the following concrete example: the `EnzoMethodMHDVlct` Method is
/// a predictor-corrector scheme that needs to compute the primitive quantities
/// (including pressure) from arrays holding integration quantities at the 
/// start of the timestep and predicted after a half timestep.
///   - Integration quantities at the start of the timestep do directly map to
///     fields of a `Block`. This means we could use technically use the old
///     version of `EnzoComputePressure` at the start of a time-step (even
///     though this would obfuscate the control-flow)
///   - However, the integration quantities at the half-timestep do not
///     directly map to the fields of a `Block`.
///   - Technically, we could have forced everything to work by making use of
///     the field-history functionallity (possibly with a minor tweak), but
///     that opens a can of worms since other Methods' expectations about the
///     number and meaning of different generations of fields are not well
///     documented. It also would have increased memory usage, obfuscated the
///     control-flow and might complicate the porting of the Method to GPUs.
/// Additionally, the old procedures would always be incompatible with computing
/// derived quantities from reconstructed fields (e.g. one might want to
/// model a spatially varying adiabatic index in a self-consistent manner with
/// Grackle)
///
/// The unorthodox implementation strategy was selected to try to minimize the
/// cost of using this wrapper. A particular goal was to facillitate the
/// `view` & `ptr_grackle` methods of `ArrayWrapper` & `BlockWrapper`
/// to be inlined. (This would not be possible if the methods were implemented
/// and called as virtual methods).

#ifndef ENZO_ENZO_FIELD_ACCESSOR_ADAPTOR_HPP
#define ENZO_ENZO_FIELD_ACCESSOR_ADAPTOR_HPP


namespace enzo_field_adaptor_detail {

  //----------------------------------------------------------------------

  class WrapperBase {};

  //----------------------------------------------------------------------

  class ArrayMapWrapper : public WrapperBase {
    /// Wrapper for EnzoEFltArrayMap.

  public:

    ArrayMapWrapper(const EnzoEFltArrayMap& array_map);

    inline CelloView<const enzo_float, 3> view(const std::string& name)
      const noexcept
    { return array_map_[name]; }
    
    inline const enzo_float* ptr_grackle(const std::string& name) const noexcept
    { return array_map_.contains(name) ? array_map_[name].data() : nullptr; }

    std::array<int,3> field_strides() const noexcept;

    void cell_width(double *hx, double *hy, double *hz) const noexcept
    { ERROR("ArrayMapWrapper::cell_width", "not implemented yet"); }

    void grackle_field_grid_props(std::array<int,3>& grid_dimension,
                                  std::array<int,3>& grid_start,
                                  std::array<int,3>& grid_end) const noexcept;

    double compute_time() const noexcept
    { ERROR("ArrayMapWrapper::compute_time", "not implemented yet"); }

  private:
    const EnzoEFltArrayMap array_map_;
  };

  //----------------------------------------------------------------------

  class BlockWrapper : public WrapperBase {
    /// Wraps a Block. This always return arrays that include ghost zones.

  public:
    BlockWrapper(Block* block, int index_history);

    inline CelloView<const enzo_float, 3> view(const std::string& name)
      const noexcept
    {
      return field_.view<enzo_float>(name,ghost_choice::include,index_history_);
    }

    inline const enzo_float* ptr_grackle(const std::string& name) const noexcept
    {
      int id_field = field_.field_id(name);

      bool correct_prec;
      switch (field_.precision(id_field)){
        case precision_default:
          correct_prec = true;
          break;
        case precision_single:
          correct_prec = std::is_same<enzo_float, float>::value;
          break;
        case precision_double:
          correct_prec = std::is_same<enzo_float, double>::value;
          break;
        case precision_quadruple:
          correct_prec = std::is_same<enzo_float, long double>::value;
          break;
        default:
          correct_prec = false;
      }

      if ((id_field >= 0) & !correct_prec){
        ERROR1("BlockWrapper::ptr_grackle",
               "%s doesn't have default precision", name.c_str());
      }
      return (const enzo_float*)field_.values(id_field, index_history_);
    }

    std::array<int,3> field_strides() const noexcept;

    void cell_width(double *hx, double *hy, double *hz) const noexcept
    { block_->cell_width(hx, hy, hz); }

    /// grid_start and grid_end will include the ghost zones
    void grackle_field_grid_props(std::array<int,3>& grid_dimension,
                                  std::array<int,3>& grid_start,
                                  std::array<int,3>& grid_end) const noexcept;

    double compute_time() const noexcept;

  private:
    Block* block_;
    Field field_;
    int index_history_;
  };

}

//----------------------------------------------------------------------

class EnzoFieldAdaptor {

  /// @class    EnzoFieldAdaptor
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Adaptor to provide a single interface for
  /// field-like access from a Block or EnzoEFltArrayMap.
  ///
  /// When this wraps a `Block` instance, this always returns data that
  /// includes the ghost zone.
  ///
  /// Instances of this class are intended to have relatively short lifetimes.
  /// It's primarily intended to be used in contexts where you are certain
  /// nothing can cause the wrapped objects to be deallocated. So a typical
  /// usage pattern case is to stack allocate an instance just for the local
  /// context where it's needed (e.g. to pass to a few functions) and then let
  /// it be destroyed. Significant care must be taken if you want this object
  /// to persist.
  ///
  /// As an aside, equivalent lifetime considerations must be accounted for
  /// when using instances of the `Field` class.
  ///
  /// @note
  /// The `EnzoFieldAdaptor::ptr_for_grackle` method should ONLY be used to
  /// pass data to Grackle. Users should never directly use that output. The
  /// first element of the pointer may only make sense when interpreted with
  /// results of `EnzoFieldAdaptor::get_grackle_field_grid_props`.

  using BlockWrapper = enzo_field_adaptor_detail::BlockWrapper;
  using ArrayMapWrapper = enzo_field_adaptor_detail::ArrayMapWrapper;

public:

  EnzoFieldAdaptor() = delete;
  EnzoFieldAdaptor(const EnzoFieldAdaptor&) = delete;
  EnzoFieldAdaptor& operator=(const EnzoFieldAdaptor&) = delete;
  EnzoFieldAdaptor(EnzoFieldAdaptor&&) = delete;
  EnzoFieldAdaptor& operator=(EnzoFieldAdaptor&&) = delete;

  /// Construct from a Block object
  ///
  /// An implicit assumption is that `block->data->field.ghosts_allocated`
  /// is always `true`, even when passing `ghost_choice::exclude`
  EnzoFieldAdaptor(Block* block, int index_history)
    : holds_block_(true),
      wrapper_(nullptr)
  {
    // allocate wrapper with placement new (to avoid heap allocation)
    wrapper_ = new(&wrapper_storage_) BlockWrapper(block, index_history);
  }

  /// Construct from a EnzoEFltArrayMap object
  EnzoFieldAdaptor(const EnzoEFltArrayMap& array_map)
    : holds_block_(false),
      wrapper_(nullptr)
  {
    // allocate wrapper with placement new (to avoid heap allocation)
    wrapper_ = new(&wrapper_storage_) ArrayMapWrapper(array_map); 
  }

  ~EnzoFieldAdaptor()
  {
    // since wrappers are allocated with placement new, we need to use explicit
    // destructors
    if (holds_block_){
      reinterpret_cast<BlockWrapper*>(wrapper_)->~BlockWrapper();
    } else {
      reinterpret_cast<ArrayMapWrapper*>(wrapper_)->~ArrayMapWrapper();
    }
  }

  /// Returns whether arr has the same strides as the wrapped fields
  inline bool consistent_with_field_strides
  (const CelloView<const enzo_float, 3>& arr) const noexcept
  {
    std::array<int,3> ref_strides = (holds_block_) ?
      reinterpret_cast<BlockWrapper*>(wrapper_)->field_strides() :
      reinterpret_cast<ArrayMapWrapper*>(wrapper_)->field_strides();

    return ( (ref_strides[0] == arr.stride(0)) &
             (ref_strides[1] == arr.stride(1)) &
             (ref_strides[2] == arr.stride(2)) );
  }

  /// Return a CelloView that acts as a view of the corresponding field
  ///
  /// When a Field is being wrapped, the ghost zones are always excluded
  /// (if we want to make this configurable in the future, the choice should be
  /// passed to the appropriate constructor)
  inline CelloView<const enzo_float, 3> view(const std::string& name)
    const noexcept
  {
    if (holds_block_){
      return reinterpret_cast<BlockWrapper*>(wrapper_)->view(name);
    } else {
      return reinterpret_cast<ArrayMapWrapper*>(wrapper_)->view(name);
    }
  }

  /// Retrieve the field as a pointer.
  ///
  /// @param name Name of the field
  /// @param require_exists Specifies behavior when `name` can't be found. When
  ///     true an error message is printed and the program aborts. Otherwise,
  ///     a `nullptr` gets returned.
  ///
  /// The returned pointers should only be used for passing data to Grackle
  ///
  /// @note
  /// You should treat the returned value as a `const enzo_float*` (do not
  /// modify the values referenced by the pointer)
  inline gr_float* ptr_for_grackle(const std::string& name,
                                   bool require_exists = false) const{
    const gr_float* ptr;
    if (holds_block_){
      ptr = reinterpret_cast<BlockWrapper*>(wrapper_)->ptr_grackle(name);
    } else {
      ptr = reinterpret_cast<ArrayMapWrapper*>(wrapper_)->ptr_grackle(name);
    }

    if ((ptr == nullptr) & (require_exists)){
      ERROR1("EnzoFieldAdaptor::ptr_for_grackle",
             "there is no array called \"%s\"", name.c_str());
    }
    return const_cast<gr_float*>(ptr);
  }

  /// Compute the grackle grid properties.
  ///
  /// This requires grid_dimension, grid_start and grid_end to be pre-allocated
  /// (they should each have space for storing 3 integers)
  inline void get_grackle_field_grid_props(int *grid_dimension,
                                           int *grid_start,
                                           int *grid_end) const noexcept
  {
    std::array<int,3> dimension_, start_, end_;

    if (holds_block_){
      reinterpret_cast<BlockWrapper*>(wrapper_)
        -> grackle_field_grid_props(dimension_, start_, end_);
    } else {
      reinterpret_cast<ArrayMapWrapper*>(wrapper_)
        -> grackle_field_grid_props(dimension_, start_, end_);
    }

    for (int i=0; i<3; i++){
      grid_dimension[i] = dimension_[i];
      grid_start[i]     = start_[i];
      grid_end[i]       = end_[i];
    }
  }

  /// Returns the cell width
  ///
  /// This is only needed for use with grackle (sometimes)
  inline void cell_width(double *hx, double *hy, double *hz) const
  {
    if (holds_block_){
      reinterpret_cast<BlockWrapper*>(wrapper_)->cell_width(hx,hy,hz);
    } else {
      reinterpret_cast<ArrayMapWrapper*>(wrapper_)->cell_width(hx,hy,hz);
    }
  }

  /// Returns the time that data was computed at
  ///
  /// This is only needed for use with grackle (for cosmological simulations)
  inline double compute_time() const
  {
    if (holds_block_){
      return reinterpret_cast<BlockWrapper*>(wrapper_)->compute_time();
    } else {
      return reinterpret_cast<ArrayMapWrapper*>(wrapper_)->compute_time();
    }
  }

private:

  /// Specifies whether the instance holds a Field or EnzoEFltArrayMap
  bool holds_block_;
  /// Pointer to wrapper around the Field or EnzoEFltArrayMap
  enzo_field_adaptor_detail::WrapperBase* wrapper_;
  /// Storage space where the wrapper_ pointer will be allocated. (This let's
  /// us avoid extra heap allocations)
  std::aligned_union<0, BlockWrapper, ArrayMapWrapper>::type wrapper_storage_;
};


#endif /* ENZO_ENZO_FIELD_ACCESSOR_ADAPTOR_HPP */
