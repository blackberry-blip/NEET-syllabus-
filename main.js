// Syllabus Data
const syllabusData = {
    physics: [
        {
            unit: "UNIT 1: PHYSICS AND MEASUREMENT",
            topics: [
                "Units of measurements, System of Units, S I Units, fundamental and derived units",
                "Least count, significant figures, Errors in measurements",
                "Dimensions of Physics quantities, dimensional analysis, and its applications"
            ]
        },
        {
            unit: "UNIT 2: KINEMATICS",
            topics: [
                "The frame of reference, motion in a straight line, Position-time graph, speed and velocity",
                "Uniform and non-uniform motion, average speed and instantaneous velocity",
                "Uniformly accelerated motion, velocity-time, position-time graph",
                "Relations for uniformly accelerated motion",
                "Scalars and Vectors, Vector Addition and subtraction, scalar and vector products",
                "Unit Vector, Resolution of a Vector, Relative Velocity",
                "Motion in a plane, Projectile Motion, Uniform Circular Motion"
            ]
        },
        {
            unit: "UNIT 3: LAWS OF MOTION",
            topics: [
                "Force and inertia, Newton's First law of motion",
                "Momentum, Newton's Second Law of motion, Impulses",
                "Newton's Third Law of motion",
                "Law of conservation of linear momentum and its applications",
                "Equilibrium of concurrent forces",
                "Static and Kinetic friction, laws of friction, rolling friction",
                "Dynamics of uniform circular motion: centripetal force and its applications"
            ]
        },
        {
            unit: "UNIT 4: WORK, ENERGY, AND POWER",
            topics: [
                "Work done by a constant force and a variable force",
                "Kinetic and potential energies, work-energy theorem, power",
                "The potential energy of spring conservation of mechanical energy",
                "Conservative and non-conservative forces",
                "Motion in a vertical circle",
                "Elastic and inelastic collisions in one and two dimensions"
            ]
        },
        {
            unit: "UNIT 5: ROTATIONAL MOTION",
            topics: [
                "Centre of the mass of a two-particle system, Centre of the mass of a rigid body",
                "Basic concepts of rotational motion; moment of a force; torque, angular momentum",
                "Conservation of angular momentum and its applications",
                "The moment of inertia, the radius of gyration",
                "Values of moments of inertia for simple geometrical objects",
                "Parallel and perpendicular axes theorems, and their applications",
                "Equilibrium of rigid bodies, rigid body rotation and equations of rotational motion",
                "Comparison of linear and rotational motions"
            ]
        },
        {
            unit: "UNIT 6: GRAVITATION",
            topics: [
                "The universal law of gravitation",
                "Acceleration due to gravity and its variation with altitude and depth",
                "Kepler's law of planetary motion",
                "Gravitational potential energy; gravitational potential",
                "Escape velocity, Motion of a satellite, orbital velocity, time period and energy of satellite"
            ]
        },
        {
            unit: "UNIT 7: PROPERTIES OF SOLIDS AND LIQUIDS",
            topics: [
                "Elastic behaviour, Stress-strain relationship, Hooke's Law",
                "Young's modulus, bulk modulus, modulus of rigidity",
                "Pressure due to a fluid column; Pascal's law and its applications",
                "Effect of gravity on fluid pressure",
                "Viscosity, Stokes' law, terminal velocity",
                "Streamline and turbulent flow, critical velocity",
                "Bernoulli's principle and its applications",
                "Surface energy and surface tension, angle of contact",
                "Excess of pressure across a curved surface, application of surface tension",
                "Heat, temperature, thermal expansion; specific heat capacity, calorimetry",
                "Change of state, latent heat",
                "Heat transfer-conduction, convection, and radiation"
            ]
        },
        {
            unit: "UNIT 8: THERMODYNAMICS",
            topics: [
                "Thermal equilibrium, zeroth law of thermodynamics, the concept of temperature",
                "Heat, work, and internal energy",
                "The first law of thermodynamics, isothermal and adiabatic processes",
                "The second law of thermodynamics: reversible and irreversible processes"
            ]
        },
        {
            unit: "UNIT 9: KINETIC THEORY OF GASES",
            topics: [
                "Equation of state of a perfect gas, work done on compressing a gas",
                "Kinetic theory of gases - assumptions, the concept of pressure",
                "Kinetic interpretation of temperature: RMS speed of gas molecules",
                "Degrees of freedom, Law of equipartition of energy",
                "Applications to specific heat capacities of gases",
                "Mean free path, Avogadro's number"
            ]
        },
        {
            unit: "UNIT 10: OSCILLATIONS AND WAVES",
            topics: [
                "Oscillations and periodic motion – time period, frequency, displacement as a function of time",
                "Periodic functions, Simple harmonic motion (S.H.M.) and its equation; phase",
                "Oscillations of a spring -restoring force and force constant",
                "Energy in S.H.M. - Kinetic and potential energies",
                "Simple pendulum - derivation of expression for its time period",
                "Wave motion, Longitudinal and transverse waves, speed of travelling wave",
                "Displacement relation for a progressive wave",
                "Principle of superposition of waves, reflection of waves",
                "Standing waves in strings and organ pipes, fundamental mode and harmonics",
                "Beats"
            ]
        },
        {
            unit: "UNIT 11: ELECTROSTATICS",
            topics: [
                "Electric charges: Conservation of charge",
                "Coulomb's law forces between two point charges, forces between multiple charges",
                "Superposition principle and continuous charge distribution",
                "Electric field: Electric field due to a point charge, Electric field lines",
                "Electric dipole, Electric field due to a dipole",
                "Torque on a dipole in a uniform electric field",
                "Electric flux, Gauss's law and its applications",
                "Electric potential and its calculation for a point charge, electric dipole and system of charges",
                "Potential difference, Equipotential surfaces",
                "Electrical potential energy of a system of two point charges and of electric dipole",
                "Conductors and insulators, Dielectrics and electric polarization",
                "Capacitors and capacitances, combination of capacitors",
                "Energy stored in a capacitor"
            ]
        },
        {
            unit: "UNIT 12: CURRENT ELECTRICITY",
            topics: [
                "Electric current, Drift velocity, mobility and their relation with electric current",
                "Ohm's law, Electrical resistance, V-I characteristics",
                "Electrical energy and power, Electrical resistivity and conductivity",
                "Series and parallel combinations of resistors",
                "Temperature dependence of resistance",
                "Internal resistance, potential difference and emf of a cell",
                "Combination of cells in series and parallel",
                "Kirchhoff's laws and their applications",
                "Wheatstone bridge, Metre Bridge"
            ]
        },
        {
            unit: "UNIT 13: MAGNETIC EFFECTS OF CURRENT AND MAGNETISM",
            topics: [
                "Biot - Savart law and its application to current carrying circular loop",
                "Ampere's law and its applications to infinitely long current carrying straight wire and solenoid",
                "Force on a moving charge in uniform magnetic and electric fields",
                "Force on a current-carrying conductor in a uniform magnetic field",
                "The force between two parallel currents carrying conductors-definition of ampere",
                "Torque experienced by a current loop in a uniform magnetic field",
                "Moving coil galvanometer, its sensitivity, and conversion to ammeter and voltmeter",
                "Current loop as a magnetic dipole and its magnetic dipole moment",
                "Bar magnet as an equivalent solenoid, magnetic field lines",
                "Magnetic field due to a magnetic dipole (bar magnet) along its axis and perpendicular to its axis",
                "Torque on a magnetic dipole in a uniform magnetic field",
                "Para-, dia- and ferromagnetic substances with examples",
                "Effect of temperature on magnetic properties"
            ]
        },
        {
            unit: "UNIT 14: ELECTROMAGNETIC INDUCTION AND ALTERNATING CURRENTS",
            topics: [
                "Electromagnetic induction: Faraday's law",
                "Induced emf and current: Lenz's Law, Eddy currents",
                "Self and mutual inductance",
                "Alternating currents, peak and RMS value of alternating current/ voltage",
                "Reactance and impedance: LCR series circuit, resonance",
                "Power in AC circuits, wattless current",
                "AC generator and transformer"
            ]
        },
        {
            unit: "UNIT 15: ELECTROMAGNETIC WAVES",
            topics: [
                "Displacement current",
                "Electromagnetic waves and their characteristics",
                "Transverse nature of electromagnetic waves",
                "Electromagnetic spectrum (radio waves, microwaves, infrared, visible, ultraviolet, X-rays, Gamma rays)",
                "Applications of e.m. waves"
            ]
        },
        {
            unit: "UNIT 16: OPTICS",
            topics: [
                "Reflection of light, spherical mirrors, mirror formula",
                "Refraction of light at plane and spherical surfaces",
                "Thin lens formula and lens maker formula",
                "Total internal reflection and its applications",
                "Magnification, Power of a Lens",
                "Combination of thin lenses in contact",
                "Refraction of light through a prism",
                "Microscope and Astronomical Telescope (reflecting and refracting) and their magnifying powers",
                "Wave optics: wavefront and Huygens' principle",
                "Laws of reflection and refraction using Huygens principle",
                "Interference, Young's double-slit experiment and expression for fringe width",
                "Coherent sources and sustained interference of light",
                "Diffraction due to a single slit, width of central maximum",
                "Polarization, plane-polarized light: Brewster's law",
                "Uses of plane-polarized light and Polaroid"
            ]
        },
        {
            unit: "UNIT 17: DUAL NATURE OF MATTER AND RADIATION",
            topics: [
                "Dual nature of radiation",
                "Photoelectric effect, Hertz and Lenard's observations",
                "Einstein's photoelectric equation: particle nature of light",
                "Matter waves-wave nature of particle, de Broglie relation"
            ]
        },
        {
            unit: "UNIT 18: ATOMS AND NUCLEI",
            topics: [
                "Alpha-particle scattering experiment; Rutherford's model of atom",
                "Bohr model, energy levels, hydrogen spectrum",
                "Composition and size of nucleus, atomic masses",
                "Mass-energy relation, mass defect",
                "Binding energy per nucleon and its variation with mass number",
                "Nuclear fission and fusion"
            ]
        },
        {
            unit: "UNIT 19: ELECTRONIC DEVICES",
            topics: [
                "Semiconductors; semiconductor diode: I-V characteristics in forward and reverse bias",
                "Diode as a rectifier; I-V characteristics of LED, the photodiode, solar cell, and Zener diode",
                "Zener diode as a voltage regulator",
                "Logic gates (OR, AND, NOT, NAND and NOR)"
            ]
        },
        {
            unit: "UNIT 20: EXPERIMENTAL SKILLS",
            topics: [
                "Vernier calipers-its use to measure the internal and external diameter and depth of a vessel",
                "Screw gauge-its use to determine thickness/ diameter of thin sheet/wire",
                "Simple Pendulum-dissipation of energy by plotting a graph between the square of amplitude and time",
                "Metre Scale - the mass of a given object by the principle of moments",
                "Young's modulus of elasticity of the material of a metallic wire",
                "Surface tension of water by capillary rise and effect of detergents",
                "Co-efficient of Viscosity of a given viscous liquid by measuring terminal velocity of a given spherical body",
                "Speed of sound in air at room temperature using a resonance tube",
                "Specific heat capacity of a given (i) solid and (ii) liquid by method of mixtures",
                "The resistivity of the material of a given wire using a metre bridge",
                "The resistance of a given wire using Ohm's law",
                "Resistance and figure of merit of a galvanometer by half deflection method",
                "The focal length of; (i) Convex mirror (ii) Concave mirror, and (iii) Convex lens, using the parallax method",
                "The plot of the angle of deviation vs angle of incidence for a triangular prism",
                "Refractive index of a glass slab using a travelling microscope",
                "Characteristic curves of a p-n junction diode in forward and reverse bias",
                "Characteristic curves of a Zener diode and finding reverse break down voltage",
                "Identification of Diode, LED, Resistor, A capacitor from a mixed collection of such items"
            ]
        }
    ],
    chemistry: [
        {
            unit: "UNIT 1: SOME BASIC CONCEPTS IN CHEMISTRY",
            topics: [
                "Matter and its nature, Dalton's atomic theory",
                "Concept of atom, molecule, element, and compound",
                "Laws of chemical combination",
                "Atomic and molecular masses, mole concept, molar mass, percentage composition",
                "Empirical and molecular formulae",
                "Chemical equations and stoichiometry"
            ]
        },
        {
            unit: "UNIT 2: ATOMIC STRUCTURE",
            topics: [
                "Nature of electromagnetic radiation, photoelectric effect",
                "Spectrum of the hydrogen atom",
                "Bohr model of a hydrogen atom - its postulates",
                "Derivation of the relations for the energy of the electron and radii of the different orbits",
                "Limitations of Bohr's model",
                "Dual nature of matter, de Broglie's relationship",
                "Heisenberg uncertainty principle",
                "Elementary ideas of quantum mechanics",
                "Quantum mechanical model of the atom, its important features",
                "Concept of atomic orbitals as one-electron wave functions",
                "Variation of Ψ and Ψ² with r for 1s and 2s orbitals",
                "Various quantum numbers (principal, angular momentum, and magnetic quantum numbers) and their significance",
                "Shapes of s, p, and d - orbitals",
                "Electron spin and spin quantum number",
                "Rules for filling electrons in orbitals – Aufbau principle, Pauli's exclusion principle and Hund's rule",
                "Electronic configuration of elements",
                "Extra stability of half-filled and completely filled orbitals"
            ]
        },
        {
            unit: "UNIT 3: CHEMICAL BONDING AND MOLECULAR STRUCTURE",
            topics: [
                "Kossel - Lewis approach to chemical bond formation",
                "Concept of ionic and covalent bonds",
                "Ionic Bonding: Formation of ionic bonds, factors affecting the formation of ionic bonds",
                "Calculation of lattice enthalpy",
                "Covalent Bonding: Concept of electronegativity, Fajan's rule, dipole moment",
                "Valence Shell Electron Pair Repulsion (VSEPR) theory and shapes of simple molecules",
                "Quantum mechanical approach to covalent bonding: Valence bond theory",
                "Concept of hybridization involving s, p, and d orbitals",
                "Resonance",
                "Molecular Orbital Theory - Its important features",
                "LCAOs, types of molecular orbitals (bonding, antibonding), sigma and pi-bonds",
                "Molecular orbital electronic configurations of homonuclear diatomic molecules",
                "Concept of bond order, bond length, and bond energy",
                "Elementary idea of metallic bonding",
                "Hydrogen bonding and its applications"
            ]
        },
        {
            unit: "UNIT 4: CHEMICAL THERMODYNAMICS",
            topics: [
                "Fundamentals of thermodynamics: System and surroundings",
                "Extensive and intensive properties, state functions, types of processes",
                "The first law of thermodynamics - Concept of work, heat internal energy and enthalpy",
                "Heat capacity, molar heat capacity",
                "Hess's law of constant heat summation",
                "Enthalpies of bond dissociation, combustion, formation, atomization, sublimation, phase transition, hydration, ionization, and solution",
                "The second law of thermodynamics - Spontaneity of processes",
                "ΔS of the universe and ΔG of the system as criteria for spontaneity",
                "ΔG° (Standard Gibbs energy change) and equilibrium constant"
            ]
        },
        {
            unit: "UNIT 5: SOLUTIONS",
            topics: [
                "Different methods for expressing the concentration of solution - molality, molarity, mole fraction, percentage (by volume and mass both)",
                "Vapour pressure of solutions and Raoult's Law",
                "Ideal and non-ideal solutions, vapour pressure - composition, plots for ideal and non-ideal solutions",
                "Colligative properties of dilute solutions - relative lowering of vapour pressure, depression of freezing point, elevation of boiling point and osmotic pressure",
                "Determination of molecular mass using colligative properties",
                "Abnormal value of molar mass, van't Hoff factor and its significance"
            ]
        },
        {
            unit: "UNIT 6: EQUILIBRIUM",
            topics: [
                "Meaning of equilibrium, the concept of dynamic equilibrium",
                "Equilibria involving physical processes: Solid-liquid, liquid - gas and solid-gas equilibria",
                "Henry's law, General characteristics of equilibrium involving physical processes",
                "Equilibrium involving chemical processes: Law of chemical equilibrium",
                "Equilibrium constants (Kp and Kc) and their significance",
                "Significance of ΔG and ΔG° in chemical equilibrium",
                "Factors affecting equilibrium concentration, pressure, temperature, effect of catalyst",
                "Le Chatelier's principle",
                "Ionic equilibrium: Weak and strong electrolytes",
                "Ionization of electrolytes, various concepts of acids and bases (Arrhenius, Bronsted - Lowry and Lewis)",
                "Acid-base equilibria (including multistage ionization) and ionization constants",
                "Ionization of water, pH scale, common ion effect",
                "Hydrolysis of salts and pH of their solutions",
                "Solubility of sparingly soluble salts and solubility products",
                "Buffer solutions"
            ]
        },
        {
            unit: "UNIT 7: REDOX REACTIONS AND ELECTROCHEMISTRY",
            topics: [
                "Electronic concepts of oxidation and reduction, redox reactions",
                "Oxidation number, rules for assigning oxidation number",
                "Balancing of redox reactions",
                "Electrolytic and metallic conduction",
                "Conductance in electrolytic solutions, molar conductivities and their variation with concentration",
                "Kohlrausch's law and its applications",
                "Electrochemical cells - Electrolytic and Galvanic cells",
                "Different types of electrodes, electrode potentials including standard electrode potential",
                "Half-cell and cell reactions, emf of a Galvanic cell and its measurement",
                "Nernst equation and its applications",
                "Relationship between cell potential and Gibbs' energy change",
                "Dry cell and lead accumulator",
                "Fuel cells"
            ]
        },
        {
            unit: "UNIT 8: CHEMICAL KINETICS",
            topics: [
                "Rate of a chemical reaction",
                "Factors affecting the rate of reactions: concentration, temperature, pressure, and catalyst",
                "Elementary and complex reactions",
                "Order and molecularity of reactions",
                "Rate law, rate constant and its units",
                "Differential and integral forms of zero and first-order reactions",
                "Characteristics and half-lives of zero and first-order reactions",
                "Effect of temperature on the rate of reactions",
                "Arrhenius theory, activation energy and its calculation",
                "Collision theory of bimolecular gaseous reactions (no derivation)"
            ]
        },
        {
            unit: "UNIT 9: CLASSIFICATION OF ELEMENTS AND PERIODICITY IN PROPERTIES",
            topics: [
                "Modern periodic law and present form of the periodic table",
                "s, p, d and f block elements",
                "Periodic trends in properties of elements atomic and ionic radii",
                "Ionization enthalpy, electron gain enthalpy, valence, oxidation states",
                "Chemical reactivity trends"
            ]
        },
        {
            unit: "UNIT 10: P-BLOCK ELEMENTS",
            topics: [
                "Group -13 to Group 18 Elements",
                "General Introduction: Electronic configuration",
                "General trends in physical and chemical properties of elements across the periods and down the groups",
                "Unique behaviour of the first element in each group"
            ]
        },
        {
            unit: "UNIT 11: d- AND f-BLOCK ELEMENTS",
            topics: [
                "Transition Elements - General introduction",
                "Electronic configuration, occurrence and characteristics",
                "General trends in properties of the first-row transition elements - physical properties",
                "Ionization enthalpy, oxidation states, atomic radii, colour",
                "Catalytic behaviour, magnetic properties, complex formation",
                "Interstitial compounds, alloy formation",
                "Preparation, properties, and uses of K₂Cr₂O₇, and KMnO₄",
                "Inner Transition Elements - Lanthanoids",
                "Electronic configuration, oxidation states, and lanthanoid contraction",
                "Actinoids - Electronic configuration and oxidation states"
            ]
        },
        {
            unit: "UNIT 12: CO-ORDINATION COMPOUNDS",
            topics: [
                "Introduction to coordination compounds",
                "Werner's theory; ligands, coordination number, denticity, chelation",
                "IUPAC nomenclature of mononuclear co-ordination compounds",
                "Isomerism in coordination compounds",
                "Bonding-Valence bond approach and basic ideas of Crystal field theory",
                "Colour and magnetic properties of coordination compounds",
                "Importance of co-ordination compounds (in qualitative analysis, extraction of metals and in biological systems)"
            ]
        },
        {
            unit: "UNIT 13: PURIFICATION AND CHARACTERISATION OF ORGANIC COMPOUNDS",
            topics: [
                "Purification - Crystallization, sublimation, distillation, differential extraction, and chromatography",
                "Principles and applications of purification methods",
                "Qualitative analysis - Detection of nitrogen, sulphur, phosphorus, and halogens",
                "Quantitative analysis (basic principles only)",
                "Estimation of carbon, hydrogen, nitrogen, halogens, sulphur, phosphorus",
                "Calculations of empirical formulae and molecular formulae",
                "Numerical problems in organic quantitative analysis"
            ]
        },
        {
            unit: "UNIT 14: SOME BASIC PRINCIPLES OF ORGANIC CHEMISTRY",
            topics: [
                "Tetravalency of carbon",
                "Shapes of simple molecules - hybridization (s and p)",
                "Classification of organic compounds based on functional groups",
                "Homologous series",
                "Isomerism - structural and stereoisomerism",
                "Nomenclature (Trivial and IUPAC)",
                "Covalent bond fission - Homolytic and heterolytic",
                "Free radicals, carbocations, and carbanions",
                "Stability of carbocations and free radicals",
                "Electrophiles and nucleophiles",
                "Electronic displacement in a covalent bond - Inductive effect, electromeric effect, resonance, and hyperconjugation",
                "Common types of organic reactions- Substitution, addition, elimination, and rearrangement"
            ]
        },
        {
            unit: "UNIT 15: HYDROCARBONS",
            topics: [
                "Classification, isomerism, IUPAC nomenclature of hydrocarbons",
                "General methods of preparation, properties, and reactions",
                "Alkanes - Conformations: Sawhorse and Newman projections (of ethane)",
                "Mechanism of halogenation of alkanes",
                "Alkenes - Geometrical isomerism",
                "Mechanism of electrophilic addition",
                "Addition of hydrogen, halogens, water, hydrogen halides (Markownikoffs and peroxide effect)",
                "Ozonolysis and polymerization",
                "Alkynes - Acidic character",
                "Addition of hydrogen, halogens, water, and hydrogen halides",
                "Polymerization of alkynes",
                "Aromatic hydrocarbons - Nomenclature",
                "Benzene - structure and aromaticity",
                "Mechanism of electrophilic substitution: halogenation, nitration",
                "Friedel - Craft's alkylation and acylation",
                "Directive influence of the functional group in mono-substituted benzene"
            ]
        },
        {
            unit: "UNIT 16: ORGANIC COMPOUNDS CONTAINING HALOGENS",
            topics: [
                "General methods of preparation, properties, and reactions",
                "Nature of C-X bond",
                "Mechanisms of substitution reactions",
                "Uses of halogen compounds",
                "Environmental effects of chloroform, iodoform freons, and DDT"
            ]
        },
        {
            unit: "UNIT 17: ORGANIC COMPOUNDS CONTAINING OXYGEN",
            topics: [
                "General methods of preparation, properties, reactions, and uses",
                "Alcohols: Identification of primary, secondary, and tertiary alcohols",
                "Mechanism of dehydration of alcohols",
                "Phenols: Acidic nature, electrophilic substitution reactions",
                "Halogenation, nitration and sulphonation of phenols",
                "Reimer - Tiemann reaction",
                "Ethers: Structure and properties",
                "Aldehyde and Ketones: Nature of carbonyl group",
                "Nucleophilic addition to >C=O group",
                "Relative reactivities of aldehydes and ketones",
                "Important reactions - Nucleophilic addition reactions (addition of HCN, NH₃, and its derivatives)",
                "Grignard reagent reactions with carbonyl compounds",
                "Oxidation and reduction (Wolf Kishner and Clemmensen)",
                "The acidity of α-hydrogen, aldol condensation",
                "Cannizzaro reaction, Haloform reaction",
                "Chemical tests to distinguish between aldehydes and Ketones",
                "Carboxylic Acids - Acidic strength and factors affecting it"
            ]
        },
        {
            unit: "UNIT 18: ORGANIC COMPOUNDS CONTAINING NITROGEN",
            topics: [
                "General methods of preparation, Properties, reactions, and uses",
                "Amines: Nomenclature, classification structure",
                "Basic character and identification of primary, secondary, and tertiary amines",
                "Diazonium Salts: Importance in synthetic organic chemistry"
            ]
        },
        {
            unit: "UNIT 19: BIOMOLECULES",
            topics: [
                "General introduction and importance of biomolecules",
                "CARBOHYDRATES - Classification",
                "Aldoses and ketoses: monosaccharides (glucose and fructose)",
                "Constituent monosaccharides of oligosaccharides (sucrose, lactose, and maltose)",
                "PROTEINS - Elementary Idea of α-amino acids",
                "Peptide bond, polypeptides",
                "Proteins: primary, secondary, tertiary, and quaternary structure (qualitative idea only)",
                "Denaturation of proteins, enzymes",
                "VITAMINS – Classification and functions",
                "NUCLEIC ACIDS – Chemical constitution of DNA and RNA",
                "Biological functions of nucleic acids",
                "Hormones (General introduction)"
            ]
        },
        {
            unit: "UNIT 20: PRINCIPLES RELATED TO PRACTICAL CHEMISTRY",
            topics: [
                "Detection of extra elements (Nitrogen, Sulphur, halogens) in organic compounds",
                "Detection of functional groups: hydroxyl (alcoholic and phenolic), carbonyl (aldehyde and ketones) carboxyl, and amino groups",
                "The chemistry involved in the preparation of inorganic compounds: Mohr's salt, potash alum",
                "The chemistry involved in the preparation of organic compounds: Acetanilide, p-nitro acetanilide, aniline yellow, iodoform",
                "The chemistry involved in the titrimetric exercises",
                "Acids, bases and the use of indicators",
                "Oxalic-acid vs KMnO₄, Mohr's salt vs KMnO₄",
                "Chemical principles involved in the qualitative salt analysis",
                "Cations – Pb²⁺, Cu²⁺, Al³⁺, Fe³⁺, Zn²⁺, Ni²⁺, Ca²⁺, Ba²⁺, Mg²⁺, NH₄⁺",
                "Anions- CO₃²⁻, S²⁻, SO₄²⁻, NO₃⁻, NO₂⁻, Cl⁻, Br⁻, I⁻ (Insoluble salts excluded)",
                "Chemical principles involved in the following experiments:",
                "1. Enthalpy of solution of CuSO₄",
                "2. Enthalpy of neutralization of strong acid and strong base",
                "3. Preparation of lyophilic and lyophobic sols",
                "4. Kinetic study of the reaction of iodide ions with hydrogen peroxide at room temperature"
            ]
        }
    ],
    biology: [
        {
            unit: "UNIT 1: DIVERSITY IN LIVING WORLD",
            topics: [
                "What is living?",
                "Biodiversity",
                "Need for classification",
                "Taxonomy & Systematics",
                "Concept of species and taxonomical hierarchy",
                "Binomial nomenclature",
                "Five kingdom classification",
                "Salient features and classification of Monera, Protista and Fungi into major groups",
                "Lichens, Viruses and Viroids",
                "Salient features and classification of plants into major groups-Algae, Bryophytes, Pteridophytes, Gymnosperms",
                "Salient features and classification of animals-nonchordate up to phyla level and chordate up to classes level"
            ]
        },
        {
            unit: "UNIT 2: STRUCTURAL ORGANISATION IN ANIMALS AND PLANTS",
            topics: [
                "Morphology and modifications",
                "Tissues",
                "Anatomy and functions of different parts of flowering plants: Root, stem, leaf, inflorescence",
                "Flower, fruit and seed",
                "Family (malvaceae, Cruciferae, leguminocae, compositea, graminae)",
                "Animal tissues",
                "Morphology, anatomy and functions of different systems (digestive, circulatory, respiratory, nervous and reproductive) of an insect (Frog)"
            ]
        },
        {
            unit: "UNIT 3: CELL STRUCTURE AND FUNCTION",
            topics: [
                "Cell theory and cell as the basic unit of life",
                "Structure of prokaryotic and eukaryotic cell",
                "Plant cell and animal cell",
                "Cell envelope, cell membrane, cell wall",
                "Cell organelles-structure and function",
                "Endomembrane system-endoplasmic reticulum, Golgi bodies, lysosomes, vacuoles",
                "Mitochondria, ribosomes, plastids, micro bodies",
                "Cytoskeleton, cilia, flagella, centrioles (ultra structure and function)",
                "Nucleus-nuclear membrane, chromatin, nucleolus",
                "Chemical constituents of living cells: Biomolecules-structure and function of proteins, carbohydrates, lipids, nucleic acids",
                "Enzymes-types, properties, enzyme action, classification and nomenclature of enzymes",
                "Cell division: Cell cycle, mitosis, meiosis and their significance"
            ]
        },
        {
            unit: "UNIT 4: PLANT PHYSIOLOGY",
            topics: [
                "Photosynthesis: Photosynthesis as a means of Autotrophic nutrition",
                "Site of photosynthesis take place",
                "Pigments involved in Photosynthesis (Elementary idea)",
                "Photochemical and biosynthetic phases of photosynthesis",
                "Cyclic and non cyclic and photophosphorylation",
                "Chemiosmotic hypothesis",
                "Photorespiration C3 and C4 pathways",
                "Factors affecting photosynthesis",
                "Respiration: Exchange gases",
                "Cellular respiration-glycolysis, fermentation (anaerobic), TCA cycle and electron transport system (aerobic)",
                "Energy relations- Number of ATP molecules generated",
                "Amphibolic pathways",
                "Respiratory quotient",
                "Plant growth and development: Seed germination",
                "Phases of Plant growth and plant growth rate",
                "Conditions of growth",
                "Differentiation, dedifferentiation and redifferentiation",
                "Sequence of developmental process in a plant cell",
                "Growth regulators- auxin, gibberellin, cytokinin, ethylene, ABA"
            ]
        },
        {
            unit: "UNIT 5: HUMAN PHYSIOLOGY",
            topics: [
                "Breathing and Respiration: Respiratory organs in animals (recall only)",
                "Respiratory system in humans",
                "Mechanism of breathing and its regulation in humans",
                "Exchange of gases, transport of gases and regulation of respiration",
                "Respiratory volumes",
                "Disorders related to respiration-Asthma, Emphysema, Occupational respiratory disorders",
                "Body fluids and circulation: Composition of blood, blood groups, coagulation of blood",
                "Composition of lymph and its function",
                "Human circulatory system-Structure of human heart and blood vessels",
                "Cardiac cycle, cardiac output, ECG, Double circulation",
                "Regulation of cardiac activity",
                "Disorders of circulatory system-Hypertension, Coronary artery disease, Angina pectoris, Heart failure",
                "Excretory products and their elimination: Modes of excretion- Ammonotelism, ureotelism, uricotelism",
                "Human excretory system-structure and fuction",
                "Urine formation, Osmoregulation",
                "Regulation of kidney function-Renin-angiotensin, Atrial Natriuretic Factor, ADH and Diabetes insipidus",
                "Role of other organs in excretion",
                "Disorders; Uraemia, Renal failure, Renal calculi, Nephritis",
                "Dialysis and artificial kidney",
                "Locomotion and Movement: Types of movement- ciliary, flagellar, muscular",
                "Skeletal muscle- contractile proteins and muscle contraction",
                "Skeletal system and its functions",
                "Joints",
                "Disorders of muscular and skeletal system-Myasthenia gravis, Tetany, Muscular dystrophy, Arthritis, Osteoporosis, Gout",
                "Neural control and coordination: Neuron and nerves",
                "Nervous system in humans-central nervous system, peripheral nervous system and visceral nervous system",
                "Generation and conduction of nerve impulse",
                "Chemical coordination and regulation: Endocrine glands and hormones",
                "Human endocrine system-Hypothalamus, Pituitary, Pineal, Thyroid, Parathyroid, Adrenal, Pancreas, Gonads",
                "Mechanism of hormone action (Elementary Idea)",
                "Role of hormones as messengers and regulators",
                "Hypo-and hyperactivity and related disorders (Common disorders e.g. Dwarfism, Acromegaly, Cretinism, goiter, exopthalmic goiter, diabetes, Addison's disease)"
            ]
        },
        {
            unit: "UNIT 6: REPRODUCTION",
            topics: [
                "Sexual reproduction in flowering plants: Flower structure",
                "Development of male and female gametophytes",
                "Pollination-types, agencies and examples",
                "Outbreeding devices",
                "Pollen-Pistil interaction",
                "Double fertilization",
                "Post fertilization events- Development of endosperm and embryo",
                "Development of seed and formation of fruit",
                "Special modes-apomixis, parthenocarpy, polyembryony",
                "Significance of seed and fruit formation",
                "Human Reproduction: Male and female reproductive systems",
                "Microscopic anatomy of testis and ovary",
                "Gametogenesis-spermatogenesis & oogenesis",
                "Menstrual cycle",
                "Fertilisation, embryo development upto blastocyst formation, implantation",
                "Pregnancy and placenta formation (Elementary idea)",
                "Parturition (Elementary idea)",
                "Lactation (Elementary idea)",
                "Reproductive health: Need for reproductive health and prevention of sexually transmitted diseases (STD)",
                "Birth control-Need and Methods, Contraception and Medical Termination of Pregnancy (MTP)",
                "Amniocentesis",
                "Infertility and assisted reproductive technologies – IVF, ZIFT, GIFT (Elementary idea for general awareness)"
            ]
        },
        {
            unit: "UNIT 7: GENETICS AND EVOLUTION",
            topics: [
                "Heredity and variation: Mendelian Inheritance",
                "Deviations from Mendelism-Incomplete dominance, Co-dominance, Multiple alleles and Inheritance of blood groups",
                "Pleiotropy",
                "Elementary idea of polygenic inheritance",
                "Chromosome theory of inheritance",
                "Chromosomes and genes",
                "Sex determination-In humans, birds, honey bee",
                "Linkage and crossing over",
                "Sex linked inheritance-Haemophilia, Colour blindness",
                "Mendelian disorders in humans-Thalassemia",
                "Chromosomal disorders in humans; Down's syndrome, Turner's and Klinefelter's syndromes",
                "Molecular basis of Inheritance: Search for genetic material and DNA as genetic material",
                "Structure of DNA and RNA",
                "DNA packaging",
                "DNA replication",
                "Central dogma",
                "Transcription, genetic code, translation",
                "Gene expression and regulation- Lac Operon",
                "Genome and human genome project",
                "DNA finger printing, protein biosynthesis",
                "Evolution: Origin of life",
                "Biological evolution and evidences for biological evolution from Paleontology, comparative anatomy, embryology and molecular evidence)",
                "Darwin's contribution, Modern Synthetic theory of Evolution",
                "Mechanism of evolution-Variation (Mutation and Recombination) and Natural Selection with examples",
                "Types of natural selection",
                "Gene flow and genetic drift",
                "Hardy-Weinberg's principle",
                "Adaptive Radiation",
                "Human evolution"
            ]
        },
        {
            unit: "UNIT 8: BIOLOGY AND HUMAN WELFARE",
            topics: [
                "Health and Disease",
                "Pathogens",
                "Parasites causing human diseases (Malaria, Filariasis, Ascariasis, Typhoid, Pneumonia, common cold, amoebiasis, ring worm, dengue, chikungunya)",
                "Basic concepts of immunology-vaccines",
                "Cancer, HIV and AIDS",
                "Adolescence, drug and alcohol abuse, Tobacco abuse",
                "Microbes in human welfare: In household food processing",
                "Industrial production, sewage treatment, energy generation",
                "Microbes as biocontrol agents and biofertilizers"
            ]
        },
        {
            unit: "UNIT 9: BIOTECHNOLOGY AND ITS APPLICATIONS",
            topics: [
                "Principles and process of Biotechnology: Genetic engineering (Recombinant DNA technology)",
                "Application of Biotechnology in health and agriculture",
                "Human insulin and vaccine production",
                "Gene therapy",
                "Genetically modified organisms-Bt crops",
                "Transgenic Animals",
                "Biosafety issues-Biopiracy and patents"
            ]
        },
        {
            unit: "UNIT 10: ECOLOGY AND ENVIRONMENT",
            topics: [
                "Organisms and environment",
                "Population interactions-mutualism, competition, predation, parasitism",
                "Population attributes-growth, birth rate and death rate, age distribution",
                "Ecosystem: Patterns, components",
                "Productivity and decomposition",
                "Energy flow",
                "Pyramids of number, biomass, energy",
                "Biodiversity and its conservation: Concept of Biodiversity",
                "Patterns of Biodiversity",
                "Importance of Biodiversity",
                "Loss of Biodiversity",
                "Biodiversity conservation",
                "Hotspots, endangered organisms, extinction",
                "Red Data Book, biosphere reserves, National parks and sanctuaries",
                "Sacred Groves"
            ]
        }
    ]
};

// DOM Elements
const loadingAnimation = document.getElementById('loadingAnimation');
const appContainer = document.getElementById('appContainer');
const heroSection = document.getElementById('heroSection');
const subjectContent = document.getElementById('subjectContent');
const backButton = document.getElementById('backButton');
const subjectTitle = document.getElementById('subjectTitle');
const unitsContainer = document.getElementById('unitsContainer');
const subjectCards = document.querySelectorAll('.subject-card');
const progressPercent = document.getElementById('progressPercent');
const completedTopicsElement = document.getElementById('completedTopics');
const totalTopicsElement = document.getElementById('totalTopics');
const streakCount = document.getElementById('streakCount');
const reminderModal = document.getElementById('reminderModal');
const startStudyingBtn = document.getElementById('startStudyingBtn');
const remindLaterBtn = document.getElementById('remindLaterBtn');

// State variables
let completedTopics = {};
let totalTopics = {};
let currentSubject = null;
let streak = 0;
let lastStudyDate = null;

// Initialize the app
function initApp() {
    // Load saved data from localStorage
    loadSavedData();
    
    // Calculate and display initial progress
    calculateProgress();
    
    // Show loading animation for 3 seconds
    setTimeout(() => {
        loadingAnimation.style.opacity = '0';
        setTimeout(() => {
            loadingAnimation.style.display = 'none';
            appContainer.style.display = 'block';
            
            // Check if user studied today
            checkStudyStatus();
        }, 500);
    }, 3000);
    
    // Set up event listeners
    setupEventListeners();
}

// Load saved data from localStorage
function loadSavedData() {
    // Load completed topics
    const savedCompletedTopics = localStorage.getItem('neetCompletedTopics');
    if (savedCompletedTopics) {
        completedTopics = JSON.parse(savedCompletedTopics);
    } else {
        // Initialize empty object for each subject
        completedTopics = {
            physics: {},
            chemistry: {},
            biology: {}
        };
    }
    
    // Load streak data
    const savedStreak = localStorage.getItem('neetStreak');
    if (savedStreak) {
        streak = parseInt(savedStreak);
        streakCount.textContent = streak;
    }
    
    // Load last study date
    const savedLastStudyDate = localStorage.getItem('neetLastStudyDate');
    if (savedLastStudyDate) {
        lastStudyDate = new Date(savedLastStudyDate);
    }
}

// Save data to localStorage
function saveData() {
    localStorage.setItem('neetCompletedTopics', JSON.stringify(completedTopics));
    localStorage.setItem('neetStreak', streak.toString());
    if (lastStudyDate) {
        localStorage.setItem('neetLastStudyDate', lastStudyDate.toISOString());
    }
}

// Calculate and display progress
function calculateProgress() {
    // Calculate total topics for each subject
    totalTopics = {
        physics: 0,
        chemistry: 0,
        biology: 0
    };
    
    let totalCompleted = 0;
    let grandTotal = 0;
    
    // Calculate for each subject
    for (const subject in syllabusData) {
        let subjectCompleted = 0;
        
        syllabusData[subject].forEach(unit => {
            unit.topics.forEach((topic, index) => {
                const topicKey = `${unit.unit}-${index}`;
                totalTopics[subject]++;
                grandTotal++;
                
                if (completedTopics[subject] && completedTopics[subject][topicKey]) {
                    subjectCompleted++;
                    totalCompleted++;
                }
            });
        });
    }
    
    // Update UI
    const progressPercentage = grandTotal > 0 ? Math.round((totalCompleted / grandTotal) * 100) : 0;
    progressPercent.textContent = progressPercentage;
    completedTopicsElement.textContent = totalCompleted;
    totalTopicsElement.textContent = grandTotal;
    
    // Update progress ring
    const circle = document.querySelector('.progress-ring-circle');
    const radius = circle.r.baseVal.value;
    const circumference = 2 * Math.PI * radius;
    const offset = circumference - (progressPercentage / 100) * circumference;
    circle.style.strokeDashoffset = offset;
}

// Check if user studied today
function checkStudyStatus() {
    const today = new Date();
    today.setHours(0, 0, 0, 0);
    
    if (!lastStudyDate) {
        // First time user - show reminder in the evening
        const now = new Date();
        if (now.getHours() >= 18) { // After 6 PM
            showReminderModal();
        } else {
            // Schedule reminder for 6 PM
            const reminderTime = new Date();
            reminderTime.setHours(18, 0, 0, 0);
            const timeUntilReminder = reminderTime - now;
            
            if (timeUntilReminder > 0) {
                setTimeout(showReminderModal, timeUntilReminder);
            }
        }
    } else {
        lastStudyDate.setHours(0, 0, 0, 0);
        
        if (lastStudyDate.getTime() === today.getTime()) {
            // User already studied today
            return;
        } else if (lastStudyDate.getTime() === today.getTime() - 86400000) {
            // User studied yesterday - continue streak
            // Show reminder in the evening if not studied yet
            const now = new Date();
            if (now.getHours() >= 18) {
                showReminderModal();
            } else {
                // Schedule reminder for 6 PM
                const reminderTime = new Date();
                reminderTime.setHours(18, 0, 0, 0);
                const timeUntilReminder = reminderTime - now;
                
                if (timeUntilReminder > 0) {
                    setTimeout(showReminderModal, timeUntilReminder);
                }
            }
        } else {
            // Broken streak
            streak = 0;
            streakCount.textContent = streak;
            saveData();
            
            // Show reminder
            showReminderModal();
        }
    }
}

// Show reminder modal
function showReminderModal() {
    reminderModal.classList.add('show');
}

// Hide reminder modal
function hideReminderModal() {
    reminderModal.classList.remove('show');
}

// Update study streak
function updateStreak() {
    const today = new Date();
    today.setHours(0, 0, 0, 0);
    
    if (!lastStudyDate) {
        // First study session
        streak = 1;
    } else {
        lastStudyDate.setHours(0, 0, 0, 0);
        const yesterday = new Date(today);
        yesterday.setDate(yesterday.getDate() - 1);
        
        if (lastStudyDate.getTime() === yesterday.getTime()) {
            // Consecutive day
            streak++;
        } else if (lastStudyDate.getTime() !== today.getTime()) {
            // Not yesterday and not today - reset streak
            streak = 1;
        }
    }
    
    lastStudyDate = new Date();
    streakCount.textContent = streak;
    saveData();
}

// Set up event listeners
function setupEventListeners() {
    // Subject card clicks
    subjectCards.forEach(card => {
        card.addEventListener('click', () => {
            currentSubject = card.dataset.subject;
            showSubjectContent(currentSubject);
        });
    });
    
    // Back button
    backButton.addEventListener('click', () => {
        heroSection.style.display = 'flex';
        subjectContent.style.display = 'none';
        calculateProgress();
    });
    
    // Reminder modal buttons
    startStudyingBtn.addEventListener('click', () => {
        updateStreak();
        hideReminderModal();
    });
    
    remindLaterBtn.addEventListener('click', () => {
        hideReminderModal();
        // Show again in 2 hours
        setTimeout(showReminderModal, 7200000);
    });
}

// Show subject content
function showSubjectContent(subject) {
    heroSection.style.display = 'none';
    subjectContent.style.display = 'block';
    
    // Set subject title
    subjectTitle.textContent = subject.charAt(0).toUpperCase() + subject.slice(1);
    
    // Clear previous content
    unitsContainer.innerHTML = '';
    
    // Add units for the selected subject
    syllabusData[subject].forEach(unit => {
        const unitCard = document.createElement('div');
        unitCard.className = 'unit-card';
        
        const unitHeader = document.createElement('div');
        unitHeader.className = 'unit-header';
        unitHeader.textContent = unit.unit;
        
        const unitContent = document.createElement('div');
        unitContent.className = 'unit-content';
        
        // Add topics
        unit.topics.forEach((topic, index) => {
            const topicKey = `${unit.unit}-${index}`;
            const isCompleted = completedTopics[subject] && completedTopics[subject][topicKey];
            
            const topicItem = document.createElement('div');
            topicItem.className = 'topic-item';
            
            const checkbox = document.createElement('input');
            checkbox.type = 'checkbox';
            checkbox.className = 'topic-checkbox';
            checkbox.checked = isCompleted;
            checkbox.addEventListener('change', () => {
                if (!completedTopics[subject]) {
                    completedTopics[subject] = {};
                }
                
                if (checkbox.checked) {
                    completedTopics[subject][topicKey] = true;
                    updateStreak();
                } else {
                    delete completedTopics[subject][topicKey];
                }
                
                topicText.classList.toggle('completed', checkbox.checked);
                saveData();
                calculateProgress();
            });
            
            const topicText = document.createElement('div');
            topicText.className = 'topic-text';
            if (isCompleted) topicText.classList.add('completed');
            topicText.textContent = topic;
            
            topicItem.appendChild(checkbox);
            topicItem.appendChild(topicText);
            unitContent.appendChild(topicItem);
        });
        
        // Toggle unit content on header click
        unitHeader.addEventListener('click', () => {
            unitContent.classList.toggle('show');
        });
        
        unitCard.appendChild(unitHeader);
        unitCard.appendChild(unitContent);
        unitsContainer.appendChild(unitCard);
    });
}

// Initialize the app when the DOM is loaded
document.addEventListener('DOMContentLoaded', initApp);
