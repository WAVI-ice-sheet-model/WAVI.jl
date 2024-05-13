struct AcousticImpedance{T <: Real, N <: Integer} 
    depth :: T
    densityGrains :: T
    densityPores :: T
    grainDiameter :: T
    porosity :: T
    KPores :: T
    KGrains :: T

    gammaP0 :: T
    gammaS0 :: T
    strainHardeningIndex :: T

    refT :: N
    refGrainDiameter :: T
    refDepth :: T
    refPorosity :: T
    frequencyHz :: N
end

"""
AcousticImpedance(; <kwargs>)


Keyword arguments
=================
- depth : depth within the subglacial sediment
- densityGrains :  density of the sediment grains
- densityPores : density of the pore water
- grainDiameter : sediment grain diameter
- porosity : sediment porosity
- KPores : bulk modulus of pore water
- KGrains : bulk modulus of mineral grains
- gammaP0 : reference compressional modulus (independent of the porosity, grain size, and depth in the sediment)
- gammaS0 : reference shear modulus (independent of the porosity, grain size, and depth in the sediment)
- strainHardeningIndex : strain hardening index
- refT : arbitrary time introduced solely to avoid awkward dimensions
- refGrainDiameter : reference sediment grain diameter
- refDepth : reference depth within the subglacial sediment
- refPorosity : reference sediment porosity
- frequencyHz : frequency
"""

function AcousticImpedance(; 
                        depth=200., 
                        densityGrains=2730., 
                        densityPores=1005., 
                        grainDiameter=0.1e-3, 
                        porosity=0.4, 
                        KPores=2.374e9, 
                        KGrains=3.6e10, 
                        gammaP0=3.888e8, 
                        gammaS0=4.588e7, 
                        strainHardeningIndex=0.0851, 
                        refT=1, refGrainDiameter=1e-3, 
                        refDepth=0.3, 
                        refPorosity=0.377, 
                        frequencyHz=100) 
                        
    return AcousticImpedance(
                            depth,
                            densityGrains,
                            densityPores,
                            grainDiameter,
                            porosity,
                            KPores,
                            KGrains,
                            gammaP0,
                            gammaS0,
                            strainHardeningIndex,
                            refT,
                            refGrainDiameter,
                            refDepth,
                            refPorosity,
                            frequencyHz)
end

"""
            updateAcousticImpedance() 

update acoustic impedance
"""
function updateAcousticImpedance!(aI,gravity,effectivePressure)
    # effectivePressure=(0:0.001:1)*1e6
    # effectivePressure=(1-aI.porosity)*(aI.densityGrains-aI.densityPores)*gravity*aI.depth

    refEffectivePressure=(1-aI.refPorosity)*(aI.densityGrains-aI.densityPores)*gravity*aI.refDepth

    gammaP=aI.gammaP0*((effectivePressure.*aI.grainDiameter)./(refEffectivePressure*aI.refGrainDiameter)).^(1/3)
    gammaS=aI.gammaS0*((effectivePressure.*aI.grainDiameter)./(refEffectivePressure*aI.refGrainDiameter)).^(2/3)

    Density0 = aI.porosity*aI.densityPores+(1-aI.porosity)*aI.densityGrains
    K0 = 1 ./(aI.porosity./aI.KPores+(1-aI.porosity)./aI.KGrains)

    c0=sqrt(K0./Density0)

    cP=c0./real.((1 .+((gammaP.+(4/3)*gammaS)./(Density0*c0^2)).*(im      *2*pi*aI.frequencyHz*aI.refT).^aI.strainHardeningIndex).^(-0.5))
    # alphaP=-(2*pi*aI.frequencyHz/c0)*imag((1 .+((gammaP.+(4/3)*gammaS)./(Density0*c0^2)).*(im      *2*pi*aI.frequencyHz*aI.refT).^aI.strainHardeningIndex).^(-0.5))

    cS=sqrt.(gammaS/Density0)*((2*pi*aI.frequencyHz*aI.refT)^(aI.strainHardeningIndex/2))/cos(aI.strainHardeningIndex*pi/4)
    # alphaS=2*pi*aI.frequencyHz*tan(aI.strainHardeningIndex*pi/4)./cS

    ZP=cP*Density0;
    ZS=cS*Density0;
    return ZP,ZS
end