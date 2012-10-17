(in-package :tc)


(defgeneric-1+caller mu1 (species temperature &key &allow-other-keys)
  (:fun-doc "Viscosity of a pure gas

SPECIES - a keyword for gas, or a lennard-jones 6/12 object
TEMPERATURE - Kelvin
REST - Other arguments for the default calculation method")
  (:default-method 'cck)
  (:documentation
   "Viscosity for a pure gas

MODEL - symbol, specifying the model to be used
SPECIES - a symbol, or an object with collision parameters
TEMPERATURE - Temperature, in Kelvin
REST - other optional parameters")
  (:sv-doc "Holds the symbol specifying the default method
   for calculating mu"))

#|(defgeneric mu% (model species temperature &key &allow-other-keys)
  (:documentation "Viscosity for a pure gas

MODEL - symbol, specifying the model to be used
SPECIES - a symbol, or an object with collision parameters
TEMPERATURE - Temperature, in Kelvin
REST - other optional parameters"))|#

(defgeneric-1+caller mu2 (species1 species2 p1 p2 temperature
			   &key &allow-other-keys)
  (:fun-doc "Viscosity of a binary gas mixture.  Depending on the used
  model, additional values may be returned by the specific method
  function.  See the documentation for the generic function mu2-1

Inputs:
SPECIES1,2 - symbols, objects with collision parameters
P1,2 - partial pressures of species1,2
TEMPERATURE - Temperature, in Kelvin
REST - other optional parameters that may be required by the default 
       method specified in *MU2-DEFAULT-METHOD*

Outputs: Mixture viscosity

")
  (:default-method 'mcc)
   (:documentation
   "Binary gas mixture viscosity

Inputs:
MODEL - symbol, specifying the model to be used (cck or mcc)
SPECIES1,2 - symbols, objects with collision parameters
P1,2 - partial pressures of species1,2
TEMPERATURE - Temperature, in Kelvin
REST - other optional parameters

Outputs: Mixture viscosity

The MCC model also returns:
mu1&2, the individual viscosities
mu1/p1 & mu2/p2, collision frequences for each species")
   (:sv-doc "Holds the symbol specifying the default metod for
   calculating mu2"))


(defgeneric-1+caller lambda-ig1 (species temperature &key &allow-other-keys)
  (:fun-doc "Thermal conductivity for a pure ideal gas

SPECIES - a symbol, or an object with collision parameters
TEMPERATURE - Temperature, in Kelvin
REST - other optional parameters")
  (:default-method 'mcc)
  (:documentation
   "Thermal conductivity for a pure gas

SPECIES - a symbol, or an object with collision parameters
TEMPERATURE - Temperature, in Kelvin
REST - other optional parameters
MODEL - symbol, specifying the model to be used")
  (:sv-doc "Holds the symbol specifying the default metod for
   calculating mu2"))

#|(defgeneric lambda% (model species temperature &key &allow-other-keys)
  (:documentation "thermal conductivity for a pure gas

MODEL - symbol, specifying the model to be used
SPECIES - a symbol, or an object with collision parameters
TEMPERATURE - Temperature, in Kelvin
REST - other optional parameters"))||#


(defgeneric-1+caller lambda-ig2 (species1 species2 p1 p2 temperature
		       &key &allow-other-keys)
  (:fun-doc "Binary gas mixture thermal conductivity

SPECIES1,2 - symbols, objects with collision parameters
P1,2 - partial pressures of species1,2
TEMPERATURE - Temperature, in Kelvin
REST - other optional parameters that may be required by the default method")
  (:default-method 'mcc)
  (:documentation "Binary gas mixture thermal conductivity

MODEL - symbol, specifying the model to be used
SPECIES1,2 - symbols, objects with collision parameters
P1,2 - partial pressures of species1,2
TEMPERATURE - Temperature, in Kelvin
REST - other optional parameters"))

(defgeneric-1+caller lambda1 (species temperature cv &key &allow-other-keys)
  (:fun-doc "Thermal conductivity of a gas with cv

SPECIES - a symbol, or an object with collision parameters
TEMPERATURE - Temperature, in Kelvin
CV - Heat capacity in J/mol
REST - other optional parameters")
  (:default-method 'mcc)
  (:documentation
   "Thermal conductivity for a pure gas

SPECIES - a symbol, or an object with collision parameters
TEMPERATURE - Temperature, in Kelvin
REST - other optional parameters
MODEL - symbol, specifying the model to be used")
  (:sv-doc "Holds the symbol specifying the default metod for
   calculating mu2"))

#|(defgeneric lambda-bgm% (model species1 species2 concentration1 temperature
			   &key &allow-other-keys)
  (:documentation "Binary gas mixture thermal conductivity

MODEL - symbol, specifying the model to be used
SPECIES1,2 - symbols, objects with collision parameters
CONCENTRATION1 - relative concentration of species1.
TEMPERATURE - Temperature, in Kelvin
REST - other optional parameters"))|#

(defgeneric-1+caller D12 (collision-params pressure temperature &key &allow-other-keys)
  (:fun-doc "Binary gas diffusion of a binary gas mixture, independent
  of concentrations

COLLISION-PARAMETERS - object with collision parameters
TEMPERATURE - Temperature, in Kelvin
PRESSURE - Pressure in Pa
REST - other optional parameters determined by MODEL")
  (:default-method 'cck)
  (:documentation "Binary gas diffusion of a binary gas mixture

MODEL - symbol, specifying the model to be used
COLLISION-PARAMETERS - object with collision parameters
TEMPERATURE - Temperature, in Kelvin
PRESSURE - Pressure in Pa
REST - other optional parameters determined by MODEL"))

(defgeneric-1+caller D12b (species1 species2  p1 p2 temperature &key &allow-other-keys)
  (:fun-doc "Binary gas diffusion of a binary gas mixture, dependent
  on concentrations

SPECIES1,2 - Potential parameters of species 1,2
TEMPERATURE - Temperature, in Kelvin
P1,2 - Partial pressures of species 1 and 2
REST - other optional parameters determined by MODEL")
  (:default-method 'mcc)
  (:documentation "Binary gas diffusion of a binary gas mixture,
  dependent on concentrations

MODEL - symbol, specifying the model to be used
SPECIES1,2 - Potential parameters of species 1,2
TEMPERATURE - Temperature, in Kelvin
P1,2 - Partial pressures of species 1 and 2 in Pa
REST - other optional parameters determined by MODEL"))

(defgeneric-1+caller Delta (species1 species2 p1 p2 temperature &key &allow-other-keys)
  (:fun-doc "Species concentration dependent correction to the binary
  diffusion coefficient D12

The binary diffusion is corrected as (/ D12 (- 1 Delta)) (Ferziger & Kapor, 7.3-40)

SPECIES1,2 - Potential parameters of species 1,2
TEMPERATURE - Temperature, in Kelvin
P1,2 - Partial pressures of species 1 and 2 in Pa
")
  (:default-method 'cck)
  (:documentation 
"Species concentration dependent correction to the binary diffusion coefficient D12

The binary diffusion is corrected as (/ D12 (- 1 Delta)) (Ferziger & Kapor, 7.3-40)

METHOD - CCK or MCC
SPECIES1,2 - Potential parameters of species 1,2.
TEMPERATURE - Temperature, in Kelvin
P1,2 - Partial pressures of species 1 and 2 in Pa

Restrictions: for the CCK model, only Lennard-Jones potentials are
accepted.  The reason is that the A*, B*, C*, E*, and F*, are now
calculated for Lennard-Jones potentials only.
"))

(defgeneric-1+caller alpha-T (potential1 potential2 p1 p2 temperature &key &allow-other-keys)
  (:documentation "Thermal diffusion coefficient of a binary gas mixture

MODEL - symbol, specifying the model to be used
POTENTIAL1,2 - Potential parameters for species 1,2
X1 - fraction of species 1
TEMPERATURE - Temperature, in Kelvin
REST - other optional parameters determined by MODEL")
  (:fun-doc "Thermal diffusion coefficient of a binary gas mixture

POTENTIAL1,2 - Potential parameters for species 1,2
X1 - fraction of species 1
TEMPERATURE - Temperature, in Kelvin
REST - other optional parameters determined by MODEL")
  (:default-method 'cck))



