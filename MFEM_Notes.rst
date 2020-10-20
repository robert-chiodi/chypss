Why are these the way they are?

Using FormSystemMatrix on a bilinear form (and providing a matrix to it to place the formed matrix) modifies the bilinear form? -> Because it calls Finalize() on it.
du_dt in (second arguments) in Mult(state_vector, du_dt) comes in uninitialized. DO NOT USE BEFORE SETTING TO VALID VALUE.

Why is const not often used? For example, in ParGridFunction::ProjectBdrCoefficient, why is Array<int>& attr not const?

When adding BoundaryIntegrator with ParLinearForm, holds ess_bdr (identifying the bdr to apply to) via pointer to the array. SO DO NOT LET IT GO OUT OF SCOPE OR MODIFY
