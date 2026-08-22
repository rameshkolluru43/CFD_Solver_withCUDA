#ifndef AMR_HPP
#define AMR_HPP

#include "definitions.h"
#include "Globals.h"

// Computes a scale-invariant gradient indicator (rho + pressure).
void AMR_Compute_Gradient_Indicator();

// Tags cells for split using AMR_Gradient_Threshold and AMR_Max_Fraction.
void AMR_Tag_Cells_For_Split();

// Splits tagged leaf quad cells into 4 children, up to maxLevels.
bool AMR_Split_Tagged_Cells(int maxLevels);

// Tags parent cells for 4->1 merge based on low child indicators.
void AMR_Tag_Cells_For_Merge();

// Merges tagged parents by deactivating the 4 children and reactivating parent.
bool AMR_Merge_Tagged_Cells();

// One AMR step: indicator -> tagging -> split/merge with hysteresis.
bool AMR_Adaptive_Step();

#endif
