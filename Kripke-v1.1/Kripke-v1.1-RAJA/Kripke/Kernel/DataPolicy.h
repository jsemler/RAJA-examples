/*
 * NOTICE
 *
 * This work was produced at the Lawrence Livermore National Laboratory (LLNL)
 * under contract no. DE-AC-52-07NA27344 (Contract 44) between the U.S.
 * Department of Energy (DOE) and Lawrence Livermore National Security, LLC
 * (LLNS) for the operation of LLNL. The rights of the Federal Government are
 * reserved under Contract 44.
 *
 * DISCLAIMER
 *
 * This work was prepared as an account of work sponsored by an agency of the
 * United States Government. Neither the United States Government nor Lawrence
 * Livermore National Security, LLC nor any of their employees, makes any
 * warranty, express or implied, or assumes any liability or responsibility
 * for the accuracy, completeness, or usefulness of any information, apparatus,
 * product, or process disclosed, or represents that its use would not infringe
 * privately-owned rights. Reference herein to any specific commercial products,
 * process, or service by trade name, trademark, manufacturer or otherwise does
 * not necessarily constitute or imply its endorsement, recommendation, or
 * favoring by the United States Government or Lawrence Livermore National
 * Security, LLC. The views and opinions of authors expressed herein do not
 * necessarily state or reflect those of the United States Government or
 * Lawrence Livermore National Security, LLC, and shall not be used for
 * advertising or product endorsement purposes.
 *
 * NOTIFICATION OF COMMERCIAL USE
 *
 * Commercialization of this product is prohibited without notifying the
 * Department of Energy (DOE) or Lawrence Livermore National Security.
 */

#ifndef KERNEL_VARIABLE_POLICY_H__
#define KERNEL_VARIABLE_POLICY_H__

#include<Kripke.h>
#include<Kripke/Directions.h>
#include<Kripke/DView.h>


/*
 * Define strongly-typed indices used in Kripke
 */
RAJA_INDEX_VALUE(IMaterial,    "IMaterial");     // Material ID
RAJA_INDEX_VALUE(ILegendre,    "ILegendre");     // Legendre expansion coefficient
RAJA_INDEX_VALUE(IMoment,      "IMoment");       // Spherical harmonic moment
RAJA_INDEX_VALUE(IDirection,   "IDirection");    // Local direction
RAJA_INDEX_VALUE(IGlobalGroup, "IGlobalGroup");  // Global energy group
RAJA_INDEX_VALUE(IGroup,       "IGroup");        // Local energy group
RAJA_INDEX_VALUE(IZone,        "IZone");         // Cannonical zone number
RAJA_INDEX_VALUE(IZoneIdx,     "IZoneIdx");      // Mapped zone index (sequential in hyperplane)
RAJA_INDEX_VALUE(IMix,         "IMix");          // Mixed element slot
RAJA_INDEX_VALUE(IZoneI,       "IZoneI");        // zone on the I boundary face
RAJA_INDEX_VALUE(IZoneJ,       "IZoneJ");        // zone on the K boundary face
RAJA_INDEX_VALUE(IZoneK,       "IZoneK");        // zone on the K boundary face



/**
 * Layout policies that don't change with nesting.
 */
struct FixedLayoutPolicy {
  using perm_ell = RAJA::PERM_JI;
  using perm_ellplus = RAJA::PERM_IJ;
  using perm_tlayout = RAJA::PERM_KJI;
  typedef DLayout<int, IDirection, IMoment> Layout_Ell;
  typedef DLayout<int, IDirection, IMoment> Layout_EllPlus;

  typedef DLayout<RAJA::Index_type, IZoneI, IZoneJ, IZoneK> TLayout_Zone;
};


/**
 * Layout policies tied directly to nesting.
 */
template<typename T>
struct NestingPolicy{};

template<>
struct NestingPolicy<NEST_DGZ_T> : public FixedLayoutPolicy {
  using perm_psi_phi =  RAJA::PERM_IJK;
  using perm_sigs = RAJA::PERM_IJKL;
  using perm_sigt = RAJA::PERM_IJ;
  using perm_face = RAJA::PERM_IJLK;
  typedef DLayout<int, IDirection, IGroup, IZone>    Layout_Psi;
  typedef DLayout<int, IMoment, IGlobalGroup, IZone> Layout_Phi;
  typedef DLayout<int, ILegendre, IGlobalGroup, IGlobalGroup, IMaterial> Layout_SigS;
  typedef DLayout<int, IGroup, IZone> Layout_SigT;
  
  typedef DLayout<int, IDirection, IGroup, IZoneJ, IZoneK> Layout_FaceI;
  typedef DLayout<int, IDirection, IGroup, IZoneI, IZoneK> Layout_FaceJ;
  typedef DLayout<int, IDirection, IGroup, IZoneI, IZoneJ> Layout_FaceK;
};

template<>
struct NestingPolicy<NEST_DZG_T> : public FixedLayoutPolicy {
  using perm_psi_phi =  RAJA::PERM_IKJ;
  using perm_sigs = RAJA::PERM_ILJK;
  using perm_sigt = RAJA::PERM_JI;
  using perm_face = RAJA::PERM_ILKJ;
  typedef DLayout<int, IDirection, IGroup, IZone>    Layout_Psi;
  typedef DLayout<int, IMoment, IGlobalGroup, IZone> Layout_Phi;
  typedef DLayout<int, ILegendre, IGlobalGroup, IGlobalGroup, IMaterial> Layout_SigS;
  typedef DLayout<int, IGroup, IZone> Layout_SigT;

  typedef DLayout<int, IDirection, IGroup, IZoneJ, IZoneK> Layout_FaceI;
  typedef DLayout<int, IDirection, IGroup, IZoneI, IZoneK> Layout_FaceJ;
  typedef DLayout<int, IDirection, IGroup, IZoneI, IZoneJ> Layout_FaceK;
};

template<>
struct NestingPolicy<NEST_GDZ_T> : public FixedLayoutPolicy {
  using perm_psi_phi =  RAJA::PERM_JIK;
  using perm_sigs = RAJA::PERM_JKIL;
  using perm_sigt = RAJA::PERM_IJ;
  using perm_face = RAJA::PERM_JILK;
  typedef DLayout<int, IDirection, IGroup, IZone>    Layout_Psi;
  typedef DLayout<int, IMoment, IGlobalGroup, IZone> Layout_Phi;
  typedef DLayout<int, ILegendre, IGlobalGroup, IGlobalGroup, IMaterial> Layout_SigS;
  typedef DLayout<int, IGroup, IZone> Layout_SigT;

  typedef DLayout<int, IDirection, IGroup, IZoneJ, IZoneK> Layout_FaceI;
  typedef DLayout<int, IDirection, IGroup, IZoneI, IZoneK> Layout_FaceJ;
  typedef DLayout<int, IDirection, IGroup, IZoneI, IZoneJ> Layout_FaceK;
};

template<>
struct NestingPolicy<NEST_GZD_T> : public FixedLayoutPolicy {
  using perm_psi_phi =  RAJA::PERM_JKI;
  using perm_sigs = RAJA::PERM_JKLI;
  using perm_sigt = RAJA::PERM_IJ;
  using perm_face = RAJA::PERM_JLKI;
  typedef DLayout<int, IDirection, IGroup, IZone>    Layout_Psi;
  typedef DLayout<int, IMoment, IGlobalGroup, IZone> Layout_Phi;
  typedef DLayout<int, ILegendre, IGlobalGroup, IGlobalGroup, IMaterial> Layout_SigS;
  typedef DLayout<int, IGroup, IZone> Layout_SigT;

  typedef DLayout<int, IDirection, IGroup, IZoneJ, IZoneK> Layout_FaceI;
  typedef DLayout<int, IDirection, IGroup, IZoneI, IZoneK> Layout_FaceJ;
  typedef DLayout<int, IDirection, IGroup, IZoneI, IZoneJ> Layout_FaceK;
};

template<>
struct NestingPolicy<NEST_ZDG_T> : public FixedLayoutPolicy {
  using perm_psi_phi = RAJA::PERM_KIJ;
  using perm_sigs = RAJA::PERM_LIJK;
  using perm_sigt = RAJA::PERM_JI;
  using perm_face = RAJA::PERM_LKIJ;

  typedef DLayout<int, IDirection, IGroup, IZone>    Layout_Psi;
  typedef DLayout<int, IMoment, IGlobalGroup, IZone> Layout_Phi;
  typedef DLayout<int, ILegendre, IGlobalGroup, IGlobalGroup, IMaterial> Layout_SigS;
  typedef DLayout<int, IGroup, IZone> Layout_SigT;

  typedef DLayout<int, IDirection, IGroup, IZoneJ, IZoneK> Layout_FaceI;
  typedef DLayout<int, IDirection, IGroup, IZoneI, IZoneK> Layout_FaceJ;
  typedef DLayout<int, IDirection, IGroup, IZoneI, IZoneJ> Layout_FaceK;
};

template<>
struct NestingPolicy<NEST_ZGD_T> : public FixedLayoutPolicy {
  using perm_psi_phi = RAJA::PERM_KJI;
  using perm_sigs = RAJA::PERM_LJKI;
  using perm_sigt = RAJA::PERM_JI;
  using perm_face = RAJA::PERM_LKJI;
  typedef DLayout<int, IDirection, IGroup,       IZone>    Layout_Psi;
  typedef DLayout<int, IMoment,    IGlobalGroup, IZone> Layout_Phi;
  typedef DLayout<int, ILegendre,  IGlobalGroup, IGlobalGroup, IMaterial> Layout_SigS;
  typedef DLayout<int, IGroup,     IZone> Layout_SigT;

  typedef DLayout<int, IDirection, IGroup,       IZoneJ,       IZoneK> Layout_FaceI;
  typedef DLayout<int, IDirection, IGroup,       IZoneI,       IZoneK> Layout_FaceJ;
  typedef DLayout<int, IDirection, IGroup,       IZoneI,       IZoneJ> Layout_FaceK;
};


/**
 * Views that have fixed policies
 */
struct FixedViewPolicy {
  typedef DView<double,     DLayout<int, IZoneI> > View_dx;
  typedef DView<double,     DLayout<int, IZoneJ> > View_dy;
  typedef DView<double,     DLayout<int, IZoneK> > View_dz;
  typedef DView<Directions, DLayout<int, IDirection> > View_Directions;
  typedef DView<double,     DLayout<int, IZone> > View_Volume;
  
  typedef DView<IZoneI,     DLayout<int, IZoneIdx> > View_IdxToI;
  typedef DView<IZoneJ,     DLayout<int, IZoneIdx> > View_IdxToJ;
  typedef DView<IZoneK,     DLayout<int, IZoneIdx> > View_IdxToK;

  typedef DView<IZone,      DLayout<int, IMix> > View_MixedToZones;
  typedef DView<IMaterial,  DLayout<int, IMix> > View_MixedToMaterial;
  typedef DView<double,     DLayout<int, IMix> > View_MixedToFraction;
  typedef DView<ILegendre,  DLayout<int, IMoment> > View_MomentToCoeff;
};

/**
 * Views with policies that vary between nestings.
 */
template<typename T>
struct ViewPolicy : public FixedViewPolicy {
  // Discrete and Moment Unknowns
  typedef DView<double, typename T::Layout_Psi> View_Psi;
  typedef DView<double, typename T::Layout_Phi> View_Phi;

  // Spatial domain face indices
  typedef DView<double, typename T::Layout_FaceI> View_FaceI;
  typedef DView<double, typename T::Layout_FaceJ> View_FaceJ;
  typedef DView<double, typename T::Layout_FaceK> View_FaceK;

  // L and L+ matrices
  typedef DView<double, typename T::Layout_Ell> View_Ell;
  typedef DView<double, typename T::Layout_EllPlus> View_EllPlus;

  // Data tables
  typedef DView<double, typename T::Layout_SigS> View_SigS;
  typedef DView<double, typename T::Layout_SigT> View_SigT;
  
#ifdef RAJA_ENABLE_OPENMP
  typedef RAJA::omp_reduce reduce_policy;
#else
  typedef RAJA::seq_reduce reduce_policy;
#endif
};


/**
 * Combined Policies for Layouts, Views.
 *
 * A convenience class: makes it easier to include in application.
 */
struct FixedDataPolicy {
  static const int memory_alignment = 64;
};

template<typename T>
struct DataPolicy : public FixedDataPolicy, public NestingPolicy<T>, public ViewPolicy<NestingPolicy<T> >
{
};

#endif
