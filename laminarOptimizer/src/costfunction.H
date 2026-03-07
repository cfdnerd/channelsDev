Info<<"cost function"<<endl;
MeanT=fvc::domainIntegrate(T).value()/area;
V=fvc::domainIntegrate(xh).value()/area-(voluse+0*Foam::max(0,0.5-opt*0.015));  

PowerDissActual=0;
for (i=0; i<nObjPatch; i++)
{
     PowerDissActual =PowerDissActual - sum(
     phi.boundaryField()[conPatchList[i]]*(p.boundaryField()[conPatchList[i]]
     + 0.5*magSqr(U.boundaryField()[conPatchList[i]]))); 
}
reduce(PowerDissActual, sumOp<scalar>());

PowerDissActualOverBase=PowerDissActual/PowerDiss0;
PowerDiss=PowerDissActualOverBase;
if(GeoDim==2)
{
 PowerDiss=PowerDiss/std::pow(mesh.V()[0],1.0/3);
}
const scalar powerDissLogged = PowerDiss*PowerDiss0*std::pow(mesh.V()[0],1.0/3);

if (Pstream::master())
{ 
   std::ofstream historyFile(optimizationHistoryFile, std::ios::app);
   if (opt > 1 && ((opt - 1) % 25 == 0))
   {
     historyFile
       << std::left
       << std::setw(8) << "Iter"
       << std::setw(16) << "PDAct"
       << std::setw(16) << "PDNrm"
       << std::setw(16) << "PDCon"
       << std::setw(16) << "V"
       << std::setw(16) << "MeanT"
       << std::setw(16) << "MaxT"
       << '\n';
   }
   historyFile
     << std::right
     << std::setw(8) << opt
     << std::scientific << std::setprecision(8)
     << std::setw(16) << powerDissLogged
     << std::setw(16) << PowerDissActualOverBase
     << std::setw(16) << PowerDiss
     << std::setw(16) << V
     << std::setw(16) << MeanT
     << std::setw(16) << maxT
     << '\n';
}
