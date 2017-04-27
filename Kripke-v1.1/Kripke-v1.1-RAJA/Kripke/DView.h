#ifndef __DOMAIN_TVIEW_H__
#define __DOMAIN_TVIEW_H__

#include<string>
#include<Kripke.h>

#include "RAJA/RAJA.hxx"

template<typename IdxLin, typename ...Idxs>
struct DLayout : public RAJA::TypedLayout<IdxLin, Idxs...>{
  using Base = RAJA::TypedLayout<IdxLin, Idxs...>;
  using IndexLinear = typename Base::IndexLinear;

  constexpr static size_t n_dims = sizeof...(Idxs);

  
  inline DLayout(Grid_Data &domain, int sdom_id) :
    Base{domain.indexSize<Idxs>(sdom_id)...}
  {}

  template<typename IS>
  inline DLayout(Grid_Data &domain, int sdom_id, IS perm) :
    Base{RAJA::make_permuted_layout<n_dims>(
            {domain.indexSize<Idxs>(sdom_id)...},
            perm)}
  {}

  template<typename ... ARGS>
  RAJA_HOST_DEVICE
  inline DLayout(ARGS ... args) :
    Base{args...}
  {}

};


template<typename DataType, typename L>
struct DView {};

template<typename DataType, typename IdxLin, typename ... Idxs>
struct DView<DataType, DLayout<IdxLin, Idxs...>> :
    public RAJA::TypedView<DataType, RAJA::Layout<sizeof...(Idxs), IdxLin>, Idxs...> {
  using Base = RAJA::TypedView<DataType, RAJA::Layout<sizeof...(Idxs), IdxLin>, Idxs...>;
  using Base::index;

  template<typename IS>
  inline DView(Grid_Data &domain, int sdom_id, DataType *ptr, IS perm) :
    Base(
        ptr,
        RAJA::make_permuted_layout<sizeof...(Idxs), IdxLin>(
            {IdxLin{domain.indexSize<Idxs>(sdom_id)}...},
            perm))
  {}

  RAJA_HOST_DEVICE RAJA_INLINE 
      IdxLin index(Idxs... args) const
  {
    return this->Base::index(args...);
  }
};
template<typename IdxLin, typename...Idxs>
constexpr size_t DLayout<IdxLin, Idxs...>::n_dims;

#if 0

template<typename POL, typename IdxI, typename R, typename BODY>
RAJA_INLINE
void dForallN_expanded(Grid_Data &domain, int sdom_id, BODY const &body, R (BODY::*mf)(IdxI) const){

  RAJA::RangeSegment seg_i = domain.indexRange<IdxI>(sdom_id);

  // Call underlying forall, extracting ranges from domain
  RAJA::forallN<POL, IdxI>(seg_i, body);
}

template<typename POL, typename IdxI, typename IdxJ, typename R, typename BODY>
RAJA_INLINE
void dForallN_expanded(Grid_Data &domain, int sdom_id, BODY const &body, R (BODY::*mf)(IdxI, IdxJ) const){

  RAJA::RangeSegment seg_i = domain.indexRange<IdxI>(sdom_id);
  RAJA::RangeSegment seg_j = domain.indexRange<IdxJ>(sdom_id);

  // Call underlying forall, extracting ranges from domain
  RAJA::forallN<POL, IdxI, IdxJ>(seg_i, seg_j, body);
}



template<typename POL, typename IdxI, typename IdxJ, typename IdxK, typename R, typename BODY>
RAJA_INLINE
void dForallN_expanded(Grid_Data &domain, int sdom_id, BODY const &body, R (BODY::*mf)(IdxI, IdxJ, IdxK) const){

  RAJA::RangeSegment seg_i = domain.indexRange<IdxI>(sdom_id);
  RAJA::RangeSegment seg_j = domain.indexRange<IdxJ>(sdom_id);
  RAJA::RangeSegment seg_k = domain.indexRange<IdxK>(sdom_id);

  // Call underlying forall, extracting ranges from domain
  RAJA::forallN<POL, IdxI, IdxJ, IdxK>(seg_i, seg_j, seg_k, body);
}


template<typename POL, typename IdxI, typename IdxJ, typename IdxK, typename IdxL, typename R, typename BODY>
RAJA_INLINE
void dForallN_expanded(Grid_Data &domain, int sdom_id, BODY const &body, R (BODY::*mf)(IdxI, IdxJ, IdxK, IdxL) const){

  RAJA::RangeSegment seg_i = domain.indexRange<IdxI>(sdom_id);
  RAJA::RangeSegment seg_j = domain.indexRange<IdxJ>(sdom_id);
  RAJA::RangeSegment seg_k = domain.indexRange<IdxK>(sdom_id);
  RAJA::RangeSegment seg_l = domain.indexRange<IdxL>(sdom_id);

  // Call underlying forall, extracting ranges from domain
  RAJA::forallN<POL, IdxI, IdxJ, IdxK, IdxL>(seg_i, seg_j, seg_k, seg_l, body);
}


template<typename POLICY, typename BODY>
RAJA_INLINE 
void dForallN(Grid_Data &domain, int sdom_id, BODY body){
  dForallN_expanded<POLICY>(domain, sdom_id, body, &BODY::operator());
}

#else


#endif



template<typename POL, typename IdxI, typename BODY>
RAJA_INLINE
void dForallN(Grid_Data &domain, int sdom_id, BODY const &body){

  RAJA::RangeSegment seg_i = domain.indexRange<IdxI>(sdom_id);

  // Call underlying forall, extracting ranges from domain
  RAJA::forallN<POL, IdxI>(seg_i, body);
}

template<typename POL, typename IdxI, typename IdxJ, typename BODY>
RAJA_INLINE
void dForallN(Grid_Data &domain, int sdom_id, BODY const &body){

  RAJA::RangeSegment seg_i = domain.indexRange<IdxI>(sdom_id);
  RAJA::RangeSegment seg_j = domain.indexRange<IdxJ>(sdom_id);

  // Call underlying forall, extracting ranges from domain
  RAJA::forallN<POL, IdxI, IdxJ>(seg_i, seg_j, body);
}



template<typename POL, typename IdxI, typename IdxJ, typename IdxK, typename BODY>
RAJA_INLINE
void dForallN(Grid_Data &domain, int sdom_id, BODY const &body){

  RAJA::RangeSegment seg_i = domain.indexRange<IdxI>(sdom_id);
  RAJA::RangeSegment seg_j = domain.indexRange<IdxJ>(sdom_id);
  RAJA::RangeSegment seg_k = domain.indexRange<IdxK>(sdom_id);

  // Call underlying forall, extracting ranges from domain
  RAJA::forallN<POL, IdxI, IdxJ, IdxK>(seg_i, seg_j, seg_k, body);
}


template<typename POL, typename IdxI, typename IdxJ, typename IdxK, typename IdxL, typename BODY>
RAJA_INLINE
void dForallN(Grid_Data &domain, int sdom_id, BODY const &body){

  RAJA::RangeSegment seg_i = domain.indexRange<IdxI>(sdom_id);
  RAJA::RangeSegment seg_j = domain.indexRange<IdxJ>(sdom_id);
  RAJA::RangeSegment seg_k = domain.indexRange<IdxK>(sdom_id);
  RAJA::RangeSegment seg_l = domain.indexRange<IdxL>(sdom_id);

  // Call underlying forall, extracting ranges from domain
  RAJA::forallN<POL, IdxI, IdxJ, IdxK, IdxL>(seg_i, seg_j, seg_k, seg_l, body);
}



#endif



