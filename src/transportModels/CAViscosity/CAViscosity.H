#ifndef CAViscosity_H
#define CAViscosity_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "PtrList.H"

namespace Foam
{
namespace viscosityModels
{

class CAViscosity : public viscosityModel
{
    dictionary              CAViscosityCoeffs_;
    word                    fieldName_;

    dimensionedScalar       nu0_;      // [m^2/s]
    PtrList<dimensionedScalar> coeffs_; // nu1, nu2, ...

    bool                    clamp_;
    dimensionedScalar       clampMin_;
    dimensionedScalar       clampMax_;

protected:
    volScalarField          nu_;

    tmp<volScalarField>     calcNu() const;

public:
    TypeName("CAViscosity");

    CAViscosity
    (
        const word& name,
        const dictionary& viscosityProperties,
        const volVectorField& U,
        const surfaceScalarField& phi
    );

    virtual ~CAViscosity() = default;

    virtual tmp<volScalarField> nu() const { return nu_; }
    virtual tmp<scalarField> nu(const label patchi) const { return nu_.boundaryField()[patchi]; }

    virtual void correct();
    virtual bool read(const dictionary& viscosityProperties);
};

} // namespace viscosityModels
} // namespace Foam

#endif
