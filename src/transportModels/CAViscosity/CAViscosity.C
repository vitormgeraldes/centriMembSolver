#include "CAViscosity.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(CAViscosity, 0);
    addToRunTimeSelectionTable(viscosityModel, CAViscosity, dictionary);
}
}

using namespace Foam;

// ---------- compute ----------
tmp<volScalarField> viscosityModels::CAViscosity::calcNu() const
{
    const volScalarField& CA = U_.mesh().lookupObject<volScalarField>(fieldName_);

    tmp<volScalarField> tnu
    (
        new volScalarField
        (
            IOobject(name_, U_.time().timeName(), U_.db(), IOobject::NO_READ, IOobject::NO_WRITE),
            U_.mesh(),
            nu0_
        )
    );
    volScalarField& nuF = tnu.ref();

    nuF = nu0_;

    forAll(coeffs_, i)
    {
        const label p = i + 1;
        nuF += coeffs_[i] * pow(CA, scalar(p));
    }

    if (clamp_)
    {
        nuF = max(nuF, clampMin_);
        nuF = min(nuF, clampMax_);
    }

    return tnu;
}

// ---------- ctor ----------
viscosityModels::CAViscosity::CAViscosity
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    CAViscosityCoeffs_(viscosityProperties.optionalSubDict(typeName + "Coeffs")),
    fieldName_(CAViscosityCoeffs_.getOrDefault<word>("fieldName", "CA")),
    nu0_("nu0", dimViscosity, CAViscosityCoeffs_),
    coeffs_(0),
    clamp_(false),
    clampMin_("clampMin", dimViscosity, 0.0),
    clampMax_("clampMax", dimViscosity, VGREAT),
    nu_
    (
        IOobject(name, U_.time().timeName(), U_.db(), IOobject::NO_READ, IOobject::AUTO_WRITE),
        U_.mesh(),
        nu0_
    )
{
    // Read coefficients
    label n = 0;
    for (label i = 1; ; ++i)
    {
        const word key("nu" + Foam::name(i));
        if (!CAViscosityCoeffs_.found(key)) break;

        const dimensionedScalar nui = CAViscosityCoeffs_.get<dimensionedScalar>(key);
        coeffs_.setSize(n+1);
        coeffs_.set(n, new dimensionedScalar(nui));
        ++n;
    }

    // Clamp
    if (CAViscosityCoeffs_.found("clamp"))
    {
        const dictionary& cd = CAViscosityCoeffs_.subDict("clamp");
        if (cd.found("min")) { clampMin_ = cd.get<dimensionedScalar>("min"); clamp_ = true; }
        if (cd.found("max")) { clampMax_ = cd.get<dimensionedScalar>("max"); clamp_ = true; }
    }

    // Compute once
    nu_ = calcNu();

    Info<< "CAViscosity: field=" << fieldName_
        << ", nu0=" << nu0_
        << ", coeffs=" << coeffs_.size()
        << (clamp_ ? ", clamped" : "")
        << nl;
}

// ---------- correct ----------
void viscosityModels::CAViscosity::correct()
{
    nu_ = calcNu();
}

// ---------- read ----------
bool viscosityModels::CAViscosity::read(const dictionary& viscosityProperties)
{
    viscosityModel::read(viscosityProperties);

    CAViscosityCoeffs_ = viscosityProperties.optionalSubDict(typeName + "Coeffs");
    fieldName_ = CAViscosityCoeffs_.getOrDefault<word>("fieldName", fieldName_);
    CAViscosityCoeffs_.readEntry("nu0", nu0_);

    PtrList<dimensionedScalar> newCoeffs;
    label n = 0;
    for (label i = 1; ; ++i)
    {
        const word key("nu" + Foam::name(i));
        if (!CAViscosityCoeffs_.found(key)) break;

        const dimensionedScalar nui = CAViscosityCoeffs_.get<dimensionedScalar>(key);
        newCoeffs.setSize(n+1);
        newCoeffs.set(n, new dimensionedScalar(nui));
        ++n;
    }
    coeffs_.transfer(newCoeffs);

    clamp_ = false;
    if (CAViscosityCoeffs_.found("clamp"))
    {
        const dictionary& cd = CAViscosityCoeffs_.subDict("clamp");
        if (cd.found("min")) { clampMin_ = cd.get<dimensionedScalar>("min"); clamp_ = true; }
        if (cd.found("max")) { clampMax_ = cd.get<dimensionedScalar>("max"); clamp_ = true; }
    }

    nu_ = calcNu();

    Info<< "CAViscosity::read(): nu0=" << nu0_
        << ", coeffs=" << coeffs_.size()
        << (clamp_ ? ", clamped" : "")
        << nl;

    return true;
}
