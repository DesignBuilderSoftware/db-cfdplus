const scalar thetaInitial = 0;
const scalar deltaTheta =
    2*Foam::constant::mathematical::pi/scalar(numberOfSegments);
const scalar thetaFinal = thetaInitial + deltaTheta;
const scalar phiInitial = Foam::constant::mathematical::pi;
scalar phiFinal;
const scalar deltaPhi = deltaTheta;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
List<scalar> patchAreas(numberOfPatches, scalar(0));
List<point> patchCentroids(numberOfPatches, point::zero);
List<scalar> patchAreasConstantPhi(numberOfStripes, scalar(0));
List<point> patchCentroidsConstantPhi(numberOfStripes, point::zero);

const scalar radius = 1.0;

if (debug)
{
    Info<< "Number of stripes: " << numberOfStripes << endl;
    Info<< "Delta theta: " << deltaTheta << endl;
    Info<< "Delta phi: " << deltaPhi << endl;
}

// Computing the areas for each different Phi
for (int i = 0; i < numberOfStripes; i++)
{
    phiFinal = phiInitial - scalar(i+1.0)*deltaPhi;

if (debug)
{
    Info<< "i = " << i
        << ", phiInitial: " << phiFinal + deltaPhi
        << ", phiFinal: " << phiFinal
        << endl;
}

    // If the patch is at the poles, that is, a triangle
    if (i == 0)
    {
        point pointA(
            radius*::cos(thetaInitial)*::sin(phiInitial),
            radius*::sin(thetaInitial)*::sin(phiInitial),
            radius*::cos(phiInitial));
        point pointB(
            radius*::cos(thetaInitial)*::sin(phiFinal),
            radius*::sin(thetaInitial)*::sin(phiFinal),
            radius*::cos(phiFinal));
        point pointC(
            radius*::cos(thetaFinal)*::sin(phiFinal),
            radius*::sin(thetaFinal)*::sin(phiFinal),
            radius*::cos(phiFinal));

if (debug)
{
    Info<< "i = " << i << endl;
    Info<< "    phiInitial: " << phiInitial
        << ", phiFinal: " << phiFinal
        << endl;
    Info<< "    Point A for 1st stripe: " << pointA << endl;
    Info<< "    Point B for 1st stripe: " << pointB << endl;
    Info<< "    Point C for 1st stripe: " << pointC << endl;
}

        patchCentroidsConstantPhi[i] = (pointA+pointB+pointC)/3.0;
        patchAreasConstantPhi[i] =
            mag
            (
                0.5
              * ((pointB-pointA)^(pointC-pointA))
            );
    }
    else if (i == numberOfStripes-1)
    {
        point pointA(
            radius*::cos(thetaInitial)*::sin(phiFinal+deltaPhi),
            radius*::sin(thetaInitial)*::sin(phiFinal+deltaPhi),
            radius*::cos(phiFinal+deltaPhi));
        point pointB(
            radius*::cos(thetaFinal)*::sin(phiFinal+deltaPhi),
            radius*::sin(thetaFinal)*::sin(phiFinal+deltaPhi),
            radius*::cos(phiFinal+deltaPhi));
        point pointC(
            radius*::cos(thetaInitial)*::sin(phiFinal),
            radius*::sin(thetaInitial)*::sin(phiFinal),
            radius*::cos(phiFinal));

if (debug)
{
    Info<< "i = " << i << endl;
    Info<< "    phiInitial: " << phiFinal+deltaPhi
        << ", phiFinal: " << phiFinal
        << endl;
    Info<< "    Point A for last stripe: " << pointA << endl;
    Info<< "    Point B for last stripe: " << pointB << endl;
    Info<< "    Point C for last stripe: " << pointC << endl;
}

        patchCentroidsConstantPhi[i] = (pointA+pointB+pointC)/3.0;
        patchAreasConstantPhi[i] =
            mag
            (
                0.5
              * ((pointB-pointA)^(pointC-pointA))
            );
    }
    else // The patch is a quadrilateral
    {
        point pointA(
            radius*::cos(thetaInitial)*::sin(phiFinal+deltaPhi),
            radius*::sin(thetaInitial)*::sin(phiFinal+deltaPhi),
            radius*::cos(phiFinal+deltaPhi));
        point pointB(
            radius*::cos(thetaFinal)*::sin(phiFinal+deltaPhi),
            radius*::sin(thetaFinal)*::sin(phiFinal+deltaPhi),
            radius*::cos(phiFinal+deltaPhi));
        point pointC(
            radius*::cos(thetaInitial)*::sin(phiFinal),
            radius*::sin(thetaInitial)*::sin(phiFinal),
            radius*::cos(phiFinal));
        point pointD(
            radius*::cos(thetaFinal)*::sin(phiFinal),
            radius*::sin(thetaFinal)*::sin(phiFinal),
            radius*::cos(phiFinal));

if (debug)
{
    Info<< "i = " << i << endl;
    Info<< "    phiInitial: " << phiFinal+deltaPhi
        << ", phiFinal: " << phiFinal
        << endl;
    Info<< "    Point A : " << pointA << endl;
    Info<< "    Point B: " << pointB << endl;
    Info<< "    Point C: " << pointC << endl;
    Info<< "    Point D: " << pointD << endl;
}

        patchCentroidsConstantPhi[i] = (pointA+pointB+pointC+pointD)/4.0;
        patchAreasConstantPhi[i] = calcAreaFrom4Points
        (
            pointA,
            pointB,
            pointD,
            pointC,
            patchCentroidsConstantPhi[i]
        );
    }
}

if (debug)
{
    Info<< "Areas for constant phi: " << patchAreasConstantPhi << endl;

    for (int i = 0; i < numberOfStripes; i++)
    {
        Info<< "Stripe " << i
            << " centroid: " << patchCentroidsConstantPhi[i]
            << endl;
    }
}

for (label i = 0, j = 0, k = 0; i < numberOfPatches; i++)
{
    scalar localTheta = thetaInitial/2.0+k*deltaTheta;

    if (i < numberOfSegments*(j+1) && i >= numberOfSegments*j)
    {
        patchAreas[i] = patchAreasConstantPhi[j];
        patchCentroids[i].x() =
            patchCentroidsConstantPhi[j].x()*Foam::cos(localTheta)
          - patchCentroidsConstantPhi[j].y()*Foam::sin(localTheta);
        patchCentroids[i].y() =
            patchCentroidsConstantPhi[j].x()*Foam::sin(localTheta)
          + patchCentroidsConstantPhi[j].y()*Foam::cos(localTheta);
        patchCentroids[i].z() = patchCentroidsConstantPhi[j].z();
        k++;
    }
    else
    {
        j++;
        i--;
        k = 0;
    }
}

if (debug)
{
    for (label i = 0; i < numberOfPatches; i++)
    {
        if ((i%numberOfSegments) == 0 || i == 0)
        {
            Info<< "Stripe " << i/numberOfSegments << " areas: " << endl;
        }

        Info<< "    " << patchAreas[i] << endl;
    }

    for (label i = 0; i < numberOfPatches; i++)
    {
        if ((i%numberOfSegments) == 0 || i == 0)
        {
            Info<< "Stripe " << i/numberOfSegments << " centroids: " << endl;
        }

        Info<< "    " << patchCentroids[i] << endl;
    }
}
