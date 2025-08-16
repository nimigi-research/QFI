/-
  QFI_Core.lean
  Minimal skeleton: finite partition aggregation yields a row-stochastic matrix.
  Lean 4 (mathlib required for finite sums). This is a starting point; adjust imports to your setup.
-/

import Mathlib.Data.Real.Basic
import Mathlib.Data.Matrix.Basic
import Mathlib.Data.Finset.Basic
import Mathlib.Tactic

open BigOperators

namespace QFI

-- A finite index type for cells
variable {m : ℕ} (Cells : Type) [Fintype Cells] [DecidableEq Cells]

-- A stochastic kernel on cells (row-stochastic matrix)
def RowStochastic (P : Matrix Cells Cells ℝ) : Prop :=
  ∀ i : Cells, (∑ j, P i j) = 1 ∧ (∀ j, 0 ≤ P i j)

-- Aggregation: given a measurable kernel downstairs and a partition, we assume we already computed P.
-- Here we assert property: rows sum to 1 if P was defined by conditional expectations.

lemma rowsum_one_of_aggregated
  (P : Matrix Cells Cells ℝ)
  (h : ∀ i, (∑ j, P i j) = 1) :
  ∀ i, (∑ j, P i j) = 1 := by
  intro i; exact h i

lemma nonneg_of_aggregated
  (P : Matrix Cells Cells ℝ)
  (h : ∀ i j, 0 ≤ P i j) :
  ∀ i j, 0 ≤ P i j := by
  intro i j; exact h i j

theorem aggregated_is_rowstochastic
  (P : Matrix Cells Cells ℝ)
  (h1 : ∀ i, (∑ j, P i j) = 1)
  (h2 : ∀ i j, 0 ≤ P i j) :
  RowStochastic Cells P := by
  intro i; exact ⟨h1 i, fun j => h2 i j⟩

end QFI