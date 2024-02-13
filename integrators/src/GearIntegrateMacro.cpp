int i, j, k;
Molecule *p_mol;
#ifdef GEARTHERMO
double local_Ekin = 0.0;
#endif
#ifdef GEARBOXRESC
double scaling = 1.0;
#endif

for (i = 0; i < Nsteps; i++)
{
#ifdef GEARTHERMO
    local_Ekin = 0.0;
#endif
// prediction (with/without Ekin calculation)
#if defined(PARALLEL) && defined(GEARTHERMO)
#pragma omp parallel num_threads(thread_count) reduction(+ : local_Ekin)
    {
#elif defined(PARALLEL)
#pragma omp parallel num_threads(thread_count)
    {
#endif
        int m, l;
#ifdef PARALLEL
#pragma omp for
#endif
        for (m = 0; m < system->noMolecules; m++)
        {
            for (l = 0; l < system->molecules[m].noAtoms; l++)
            {
                *system->molecules[m].atoms[l].R = (*A) * (*system->molecules[m].atoms[l].R);
#ifdef GEARTHERMO
                local_Ekin += system->molecules[m].atoms[l].mass * system->molecules[m].atoms[l].R->RowRowDotProduct(predVel, predVel); // division by h^2 done post in AbstractGearIntegrator
#endif
            }
        }
#ifdef PARALLEL
    }
#endif
#ifdef GEARTHERMO
    local_Ekin /= h * h;
    Tinst = local_Ekin / system->noDegrOfFred;
#endif

// prediction for ExtDoFs
#ifdef GEAREXTDOFS
    *ExtDoFs = (*A) * (*ExtDoFs);
#endif

// rescaling to predicted values
#ifdef GEARBOXRESC
    scaling = GetScalingFactor(1);
    if (scaling != 1.0)
    {
        system->RescaleBox(scaling, positionsAt, molecular_based);
    }
#endif

    // forces calculation
    system->CalculateForces();

// scaling factor for SHAKE (global)
#ifdef GEARBOXRESC
    lambda_Shake = GetShakeScaling();
#endif

// temperature scaling (global for System – before RHS calculation)
#ifdef GEARTHERMO
    lambda_T = GetTempScaling();
#endif

// cycle through molecules parallelized
#ifdef PARALLEL
#pragma omp parallel for num_threads(thread_count) private(j, k, p_mol)
#endif
    for (j = 0; j < system->noMolecules; j++)
    {
        p_mol = &(system->molecules[j]);
        for (k = 0; k < p_mol->noAtoms; k++)
        {
            // RHS for physical degrees of freedom
            CalculateRHS(&(p_mol->atoms[k]));
        }
        // SHAKE if constraints present
        if (p_mol->noConstrBonds > 0)
        {
            ShakeHook(p_mol);
        }
        // correction (final correction)
        Correction(p_mol);
    }

#ifdef GEAREXTDOFS
    for (j = 0; j < ExtDoFs->GetNumberOfCols(); j++)
    {
        // calculate error of j-th ExtDOF
        GForExtDoFs->operator()(j) = RHSForExtDoF(j) - ExtDoFs->operator()(2, j);
        // correct j-th ExtDOF
        *ExtDoFs = *ExtDoFs + CalculateOuterProduct(*r, *GForExtDoFs);
        // clear G for extDOFs
        GForExtDoFs->Clear();
    }
#endif

#ifdef GEARBOXRESC
    scaling = GetScalingFactor(2);
    if (scaling != 1.0)
    {
        system->RescaleBox(scaling, positionsAt, molecular_based);
    }
#endif

    if (system->boundaryCond > 0)
    {
        system->ApplyPeriodicBC(positionsAt);
    }
}

// save ExtDoFs to system, calculate Eextra
#ifdef GEAREXTDOFS
system->Xi = ExtDoFs->operator()(0, 0);
system->Xivel = ExtDoFs->operator()(1, 0) / h;
if (ExtDoFs->GetNumberOfCols() > 1)
{
    system->Lambda = ExtDoFs->operator()(0, 1);
    system->Lambdavel = ExtDoFs->operator()(1, 1) / h;
}
system->Eextra = EnergyOfExtDoF();
#endif

return 0;