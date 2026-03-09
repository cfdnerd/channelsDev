//Author: Yu Minghao    Updated: May 2020 

static char help[] = "topology optimization of fluid problem\n";
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "wallDist.H"
#include "simpleControl.H"
#include "fvOptions.H"//
#include "MMA/MMA.h"
#include <cstring>
#include <iomanip>
#include <sstream>
#include <vector>
#include <diff.c>

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFvOptions.H"//
    #include "createFields.H"
    #include "readTransportProperties.H" 
    #include "initContinuityErrs.H"
    #include "readThermalProperties.H" 
    #include "readTurbulenceProperties.H"
    #include "opt_initialization.H"
    #include "TurbulenceCoupling.H"
    bool optimizerFailureStop(false);
    while (simple.loop(runTime))
    {
        #include "solverConvergenceReset.H"
        #include "update.H"
        #include "CustomKOmegaSSTModel.H"
        updateTurbulenceStageBState();
        #include "Primal_U.H"
        turbulenceKNonlinearChange = 0.0;
        turbulenceOmegaNonlinearChange = 0.0;
        turbulenceNutNonlinearChange = 0.0;
        turbulenceNonlinearConverged = !customTurbulenceActive;
        turbulenceCorrectorsUsed = 0;
        if (customTurbulenceActive)
        {
            for (label turbulenceCorrector = 0; turbulenceCorrector < turbulenceCorrectors; ++turbulenceCorrector)
            {
                const scalarField kPrevious(k.internalField());
                const scalarField omegaPrevious(omega.internalField());
                const scalarField nutPrevious(nut.internalField());

                #include "CustomKOmegaSSTModel.H"
                updateTurbulenceStageBState();
                #include "Primal_k.H"
                #include "CustomKOmegaSSTModel.H"
                updateTurbulenceStageBState();
                #include "Primal_omega.H"
                #include "CustomKOmegaSSTModel.H"
                updateTurbulenceStageBState();

                turbulenceKNonlinearChange =
                  computeRelativeFieldChange(k, kPrevious, Foam::max(kMin.value(), 1e-8));
                turbulenceOmegaNonlinearChange =
                  computeRelativeFieldChange(omega, omegaPrevious, Foam::max(omegaMin.value(), 1.0));
                turbulenceNutNonlinearChange =
                  computeRelativeFieldChange(nut, nutPrevious, Foam::max(nutFloor.value(), 1e-10));
                turbulenceCorrectorsUsed = turbulenceCorrector + 1;

                turbulenceNonlinearConverged =
                (
                    turbulenceCorrectorsUsed >= minTurbulenceCorrectors
                 && turbulenceKNonlinearChange <= turbulenceKNonlinearTol
                 && turbulenceOmegaNonlinearChange <= turbulenceOmegaNonlinearTol
                 && turbulenceNutNonlinearChange <= turbulenceNutNonlinearTol
                );

                if (turbulenceNonlinearConverged)
                {
                    break;
                }
            }
        }
        #include "Primal_T.H"
        const bool optimizationCycleUseTaThermalAdjoint =
            customTurbulenceActive
         && optimizationCycleTurbulentOptimizationMode
         && optThermalAdjointMode == word("ta");
        const volScalarField& optimizationCycleThermalAdjointField =
            optimizationCycleUseTaThermalAdjoint ? Ta : Tb;
        if (customTurbulenceActive && adjointCycleThermalAdjointOnly)
        {
            #include "Adjoint_Ta.H"
            #include "solverConvergenceWrite.H"
            #include "adjointCycleThermalWrite.H"
        }
        else if (customTurbulenceActive && adjointCycleFullAdjointVerification)
        {
            #include "Adjoint_Ta.H"
            #include "AdjointFlow_Ua.H"
            #include "Adjoint_k.H"
            #include "Adjoint_omega.H"
            #include "solverConvergenceWrite.H"
            #include "adjointCycleThermalWrite.H"
        }
        else if (customTurbulenceActive && optimizationCycleTurbulentOptimizationMode)
        {
            if (optimizationCycleUseTaThermalAdjoint)
            {
                #include "Adjoint_Ta.H"
            }
            else
            {
                #include "AdjointHeat_Tb.H"
            }
            #include "AdjointHeat_Ub.H"
            #include "AdjointFlow_Ua.H"
            #include "Adjoint_k.H"
            #include "Adjoint_omega.H"
            #include "costfunction.H"
            #include "sensitivity.H"
        }
        else if (customTurbulenceActive && turbulenceProModeActive)
        {
            #include "AdjointHeat_Tb.H"
            #include "AdjointHeat_Ub.H"
            #include "AdjointFlow_Ua.H"
            if (proModeClosureAdjointSolvesActive)
            {
                #include "Adjoint_k.H"
                #include "Adjoint_omega.H"
            }
            #include "costfunction.H"
            #include "sensitivity.H"
        }
        else
        {
            #include "AdjointHeat_Tb.H"
            #include "AdjointHeat_Ub.H"
            #include "AdjointFlow_Ua.H"
            #include "costfunction.H"              
            #include "sensitivity.H"
        }
        runTime.write();
        if (stopOptimization)
        {
            Info<< "Optimizer stop criteria met at iteration " << opt
                << " (" << optimizerStopReason << ")." << endl;
            optimizerFailureStop = optimizerStopReason != word("converged");
            break;
        }
    }
    if (optimizerFailureStop)
    {
        Info<< nl
            << "==============================================================================" << nl
            << "                      OPTIMIZER FAILURE SUMMARY" << nl
            << "-------------------------------------------------------------------------------" << nl
            << "  Status    : Oops! [Optimizer failed]" << nl
            << "  Iteration : " << opt << nl
            << "  Reason    : " << optimizerStopReason << nl
            << "==============================================================================" << endl;
    }
    #include "finalize.H"
    return 0;
}
