/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Changes, Authors and Copyright
    Based on OpenFOAM's old Mach utility.
    Modified and developed for blueCFD(R)-Kernel
    Copyright (C) 2016-2018 FSD blueCAPE Lda  http://www.bluecape.com.pt
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    calcPMVAndPPD

Description
    Calculates the local thermal comfort variable pvm (predicted mean vote),
    ppd (percent people dissatisfied) and the comfort index at each time. The
    comfort index classification is:
        1 - very cold, danger
        2 - cold, shivering
        3 - cool, unpleasant
        4 - cool, acceptable
        5 - slightly cool/acceptable
        6 - comfortable,pleasant/cool
        7 - comfortable,pleasant
        8 - comfortable,pleasant/warm
        9 - slightly warm/acceptable
        10 - warm, acceptable
        11 - warm, unpleasant
        12 - hot, very uncomfortable
        13 - very hot, danger
        14 - unoccupied
        15 - non-sedentary


    Also calculates the Operative Temperature:

        To = Tmr + (Ta * sqrt(10*U))
             -----------------------
               1 + sqrt(10*U)


    The -nowrite option just outputs the max value without writing the field.
    The -debug option just outputs the values of the variables used to compute
    the PMV and PPD for the first cell.

    The temperature conversion from Kelvin to Celsius is 273 instead of 273.15
    because it appears that the equations were designed for this reference.

    This utility computes the PMV, PPD and comfort index at the walls which
    although they are not needed it does due to the way that OpenFOAM operates
    for the whole cell. What matters is the internal field!

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//Reference value used on the original references of the code
const scalar KelvinCelsius = 273.00;

void ComfortIndex(volScalarField & rM_comfort, const volScalarField & rPMV)
{
    forAll(rPMV, celli)
    {
        if(rPMV[celli] < -3.5) rM_comfort[celli] =  1;
        if(rPMV[celli] >= -3.5 && rPMV[celli] < -2.5) rM_comfort[celli] =2;
        if(rPMV[celli] >= -2.5 && rPMV[celli] < -1.5) rM_comfort[celli] =3;
        if(rPMV[celli] >= -1.5 && rPMV[celli] < -1.0) rM_comfort[celli] =4;
        if(rPMV[celli] >= -1.0 && rPMV[celli] < -0.5) rM_comfort[celli] =5;
        if(rPMV[celli] >= -0.5 && rPMV[celli] < -0.2) rM_comfort[celli] =6;
        if(rPMV[celli] >= -0.2 && rPMV[celli] <  0.2) rM_comfort[celli] =7;
        if(rPMV[celli] >=  0.2 && rPMV[celli] <  0.5) rM_comfort[celli] =8;
        if(rPMV[celli] >=  0.5 && rPMV[celli] <  1.0) rM_comfort[celli] =9;
        if(rPMV[celli] >=  1.0 && rPMV[celli] <  1.5) rM_comfort[celli] =10;
        if(rPMV[celli] >=  1.5 && rPMV[celli] <  2.5) rM_comfort[celli] =11;
        if(rPMV[celli] >=  2.5 && rPMV[celli] <  3.5) rM_comfort[celli] =12;
        if(rPMV[celli] >=  3.5 && rPMV[celli] <  9.5) rM_comfort[celli] = 13;
        if(rPMV[celli] >= 9.5) rM_comfort[celli] = 13;
    }
}

scalar solveXN
(
    const scalar XNi,
    const scalar HCF,
    const scalar TAA,
    const scalar p2,
    const scalar P3,
    const scalar P4,
    const scalar P5,
    scalar & HC,
    int i,
    const bool debug
)
{
    scalar XN=XNi;     // XN = TCLA / 100
    scalar XF=XNi*2.0; // XF = TCLA / 50
    const scalar EPS=0.0015;

    int j = 0;

    if (debug)
    {
        Pout<<nl<<nl<<"celli: " << i << endl;
    }

    while(fabs(XN-XF) > EPS)
    {
        if (debug)
        {
            Pout<<"XF: " << XF << endl;
            Pout<<"XN: " << XN << endl;
        }

        XF=(XF+XN)/2.0;

        if (debug)
        {
            Pout<<"XF updated: " << XF << endl;
            Pout<<"100*XF- TAA: " << 100*XF- TAA << endl;
        }

        scalar HCN=2.38*Foam::pow(fabs(100*XF- TAA), 0.25);

        HC=(HCF > HCN) ? (HCF) : (HCN);

        if (debug)
        {
            Pout<<"100+P3*HC: " << 100+P3*HC << endl;
        }

        XN =
            (
                P5
              + P4*HC
              - p2*(Foam::pow(XF, 4))
            )
          / (100+P3*HC);

        if (debug)
        {
            Pout<<"XN updated: " << XN << endl;
        }
        j++;
    }

    if (debug)
    {
        Pout<<"XF: " << XF << endl;
        Pout<<"XN: " << XN << endl;
        Pout<<"HC: " << HC << endl;
        Pout<<"HCF: " << HCF << endl;
    }

    return XN;
}

void calcThermalComfortParameters
(
    const scalar & clothingIndex,
    const scalar & activityLevel,
    tmp<volScalarField> & relativityHumidity,
    const volVectorField & velocityField,
    const volScalarField & temperatureFieldKelvin,
    const volScalarField & temperatureMRTFieldKelvin,
    volScalarField & PMVField,
    volScalarField & PPDField,
    volScalarField & rComfortIndex,
    const bool debug
)
{
    const label dbgIndex=0;

    dimensionedScalar invT("invT", dimTemperature, scalar(1.0));
    const volScalarField TAA(temperatureFieldKelvin/invT);
    const volScalarField TRA(temperatureMRTFieldKelvin/invT);
    dimensionedScalar invU("invU", dimVelocity, scalar(1.0));

    volScalarField FNPS(exp(16.6536-4030.183/(TAA-KelvinCelsius+235)));

    volScalarField PA(relativityHumidity*10*FNPS);
    scalar ICL = 0.155*clothingIndex;

    scalar m = activityLevel*58.15;

    scalar FCL = (ICL < 0.078) ? (1+1.29*ICL) : (1.05+0.645*ICL);

    volScalarField HCF(12.1*Foam::sqrt(mag(velocityField)/invU));

    volScalarField TCLA
    (
        TAA
      + (35.5+KelvinCelsius-TAA)/(3.5*(6.45*ICL+0.1))
    );

    scalar p1 = ICL*FCL;
    scalar p2 = p1*3.96;
    scalar P3 = p1*100;
    volScalarField P4(p1*TAA);
    volScalarField P5(308.7-0.028*m+p2*Foam::pow((TRA/100.0), 4));

    volScalarField XN(TCLA/100.0);
    volScalarField HC(0.0*FNPS);

    if(debug)
    {
        Pout<<"dbgIndex: " << dbgIndex << endl;
        Pout<<"TAA[dbgIndex]-KelvinCelsius: " << TAA[dbgIndex]-KelvinCelsius << endl;
        Pout<<"FNPS[dbgIndex]: "           << FNPS[dbgIndex] << endl;
        Pout<<"PA[dbgIndex]: "             << PA[dbgIndex] << endl;
        Pout<<"clothing index: " << clothingIndex << endl;
        Pout<<"ICL: "            << ICL << endl;
        Pout<<"m: "              << m << endl;
        Pout<<"FCL: "            << FCL << endl;
        Pout<<"TAA[dbgIndex]: "         << TAA[dbgIndex] << endl;
        Pout<<"TRA[dbgIndex]: "         << TRA[dbgIndex] << endl;
        Pout<<"TCLA[dbgIndex]: "        << TCLA[dbgIndex] << endl;
        Pout<<"p1: " << p1 << endl;
        Pout<<"p2: " << p2 << endl;
        Pout<<"P3: " << P3 << endl;
        Pout<<"P4: " << P4[dbgIndex] << endl;
        Pout<<"P5: " << P5[dbgIndex] << endl;
    }

    forAll(XN, celli)
    {
        XN[celli] = solveXN
        (
            XN[celli],
            HCF[celli],
            TAA[celli],
            p2,
            P3,
            P4[celli],
            P5[celli],
            HC[celli],
            celli,
            debug && celli==dbgIndex
        );
    }

    volScalarField TCL(100*XN-273);

    // Skin diff loss
    volScalarField HL1(3.05*0.001*(5733-6.99*m-PA));

    // Sweat loss
    scalar HL2 = (m > 58.15) ? (0.42*(m-58.15)) : (0.0);

    // Latent respiration loss
    volScalarField HL3(1.7*0.00001*m*(5867-PA));

    //Dry respiration loss
    volScalarField HL4(0.0014*m*(34-TAA+KelvinCelsius));

    // Radiation loss
    volScalarField HL5(3.96*FCL*(pow(XN, 4)-pow((TRA/100), 4)));

    // Convection loss
    volScalarField HL6(FCL*HC*(TCL-TAA+KelvinCelsius));

    // Thermal sensation to skin tran coef
    scalar Ts = 0.303*Foam::exp(-0.036*m)+0.028;

    PMVField = Ts*(m - HL1 - HL2 - HL3 - HL4 - HL5 - HL6);

    // calculate PPD
    PPDField =
        100 - 95*Foam::exp(-0.03353*pow(PMVField, 4)-0.2179*pow(PMVField, 2));

    ComfortIndex(rComfortIndex, PMVField);

    if(debug)
    {
        Pout<<"TCL[dbgIndex]: "  << TCL[dbgIndex] << endl;
        Pout<<"HL1[dbgIndex]: "  << HL1[dbgIndex] << endl;
        Pout<<"HL2: "            << HL2 << endl;
        Pout<<"HL3[dbgIndex]: "  << HL3[dbgIndex] << endl;
        Pout<<"HL4[dbgIndex]: "  << HL4[dbgIndex] << endl;
        Pout<<"HL5[dbgIndex]: "  << HL5[dbgIndex] << endl;
        Pout<<"    FCL: "        << FCL << endl;
        Pout<<"    XN[dbgIndex]: "  << XN[dbgIndex] << endl;
        Pout<<"    TRA[dbgIndex]: " << TRA[dbgIndex] << endl;
        Pout<<"HL6[dbgIndex]: "     << HL6[dbgIndex] << endl;
        Pout<<"    FCL: "           << FCL << endl;
        Pout<<"    HC[dbgIndex]: "  << HC[dbgIndex] << endl;
        Pout<<"    TCL[dbgIndex]: " << TCL[dbgIndex] << endl;
        Pout<<"    TAA[dbgIndex]: " << TAA[dbgIndex] << endl;
        Pout<<"Ts: "                << Ts << endl;
        Pout<<"PMVField[dbgIndex]: " << PMVField[dbgIndex] << endl;
        Pout<<"PPDField[dbgIndex]: " << PPDField[dbgIndex] << endl;
        Pout<<"rComfortIndex[dbgIndex]: " << rComfortIndex[dbgIndex] << endl;
    }
}

Foam::tmp<Foam::volScalarField> retrieveHumidityField
(
    const Foam::fvMesh& mesh,
    const Foam::dictionary& PMVAndPPDDict
)
{
    if(PMVAndPPDDict.found("relativeHumidity"))
    {
        const dimensionedScalar uniformRelativeHumidity
        ("relativeHumidity", dimless, PMVAndPPDDict);

        const word humidityField("uniformRelativeHumidity");

        Info<< "Creating field " << humidityField << endl;
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    humidityField,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                uniformRelativeHumidity
            )
        );
    }
    else
    {
        const word relHumidityName("relativeHumidity");

        Info<< "Reading field " << relHumidityName << endl;
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    relHumidityName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            )
        );
    }
}

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"
    #include "addDictOption.H"

    argList::addBoolOption("debug");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    const bool debug = args.optionFound("debug");
    const word dictName("calcPMVAndPPDDict");

    // Read meshing dictionary
    IOdictionary PMVAndPPDDict
    (
       IOobject
       (
            "calcPMVAndPPDDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
       )
    );

    const scalar clothingIndex = PMVAndPPDDict.lookup<scalar>("clothingIndex");
    const scalar activityLevel = PMVAndPPDDict.lookup<scalar>("activityLevel");
    tmp<volScalarField> relativityHumidity =
        retrieveHumidityField(mesh, PMVAndPPDDict);

    instantList timeDirs = timeSelector::select0(runTime, args);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        runTime.functionObjects().start();

        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        #include "createFields.H"

        calcThermalComfortParameters(
            clothingIndex,
            activityLevel,
            relativityHumidity,
            U,
            T,
            MRT,
            PMV,
            PPD,
            comfortIndex,
            debug
            );

        // Calculate the Operative Temperature in Celsius
        {
            dimensionedScalar invU("invU", dimVelocity, scalar(1.0));
            dimensionedScalar invT("invT", dimTemperature, scalar(1.0));
            OpTempCelsius = invT *
            (
                (MRT/invT-KelvinCelsius) + ((T/invT-KelvinCelsius)*sqrt(10*mag(U/invU)))
            )
            /
            (
                scalar(1.0) + sqrt(10*mag(U/invU))
            );
        }

        Info<< endl;

        if(debug)
        {
            Pout<< "\nThermal comfort paramaters:" << endl;
            Pout<< "                  Min        Max" << endl;
            Pout<< " PMV             " << min(PMV)
                << "        " << max(PMV)
                << endl;
            Pout<< " PPD:            " << min(PPD)
                << "        " << max(PPD)
                << endl;
            Pout<< " Comfort Index:  " << min(comfortIndex)
                << "        " << max(comfortIndex)
                << endl;
        }

        PMV.write();
        PPD.write();
        comfortIndex.write();
        OpTempCelsius.write();

        runTime.functionObjects().execute();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s" << endl;
    }

    runTime.functionObjects().end();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "End\n" << endl;
}

// ************************************************************************* //
