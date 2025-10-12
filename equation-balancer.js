class ChemicalEquationBalancer {
    constructor() {
        this.elements = [];
        this.compounds = [];
        this.matrix = [];
        this.solution = [];
    }

    // Parse chemical formula to extract elements and their counts
    parseFormula(formula) {
        // Convert subscript numbers to regular numbers
        const normalizedFormula = formula
            .replace(/₀/g, '0').replace(/₁/g, '1').replace(/₂/g, '2')
            .replace(/₃/g, '3').replace(/₄/g, '4').replace(/₅/g, '5')
            .replace(/₆/g, '6').replace(/₇/g, '7').replace(/₈/g, '8')
            .replace(/₉/g, '9');

        return this.parseFormulaWithParentheses(normalizedFormula);
    }

    parseFormulaWithParentheses(formula) {
        const elements = {};
        const stack = [{}]; // Stack to handle nested parentheses
        let i = 0;

        while (i < formula.length) {
            const char = formula[i];

            if (char === '(') {
                // Start a new group
                stack.push({});
                i++;
            } else if (char === ')') {
                // End current group and apply multiplier
                i++;
                let multiplier = '';

                // Read the number after closing parenthesis
                while (i < formula.length && /\d/.test(formula[i])) {
                    multiplier += formula[i];
                    i++;
                }

                const mult = multiplier ? parseInt(multiplier) : 1;
                const group = stack.pop();
                const currentLevel = stack[stack.length - 1];

                // Add elements from group to current level with multiplier
                for (const [element, count] of Object.entries(group)) {
                    currentLevel[element] = (currentLevel[element] || 0) + (count * mult);
                }
            } else if (/[A-Z]/.test(char)) {
                // Parse element
                let element = char;
                i++;

                // Check for lowercase letter (like Cl, Ca, etc.)
                if (i < formula.length && /[a-z]/.test(formula[i])) {
                    element += formula[i];
                    i++;
                }

                // Parse number after element
                let number = '';
                while (i < formula.length && /\d/.test(formula[i])) {
                    number += formula[i];
                    i++;
                }

                const count = number ? parseInt(number) : 1;
                const currentLevel = stack[stack.length - 1];
                currentLevel[element] = (currentLevel[element] || 0) + count;
            } else {
                // Skip unknown characters
                i++;
            }
        }

        return stack[0];
    }

    // Parse the entire equation
    parseEquation(equation) {
        // Clean and normalize the equation
        equation = equation.replace(/\s+/g, ' ').trim();

        // Support different arrow notations
        const arrowRegex = /(-?>|=>?|=)/;
        if (!arrowRegex.test(equation)) {
            throw new Error('Invalid equation format. Please use ->, =>, or = as arrow.');
        }

        const parts = equation.split(arrowRegex);
        if (parts.length !== 3) {
            throw new Error('Invalid equation format. Expected format: reactants -> products');
        }

        const reactants = parts[0].trim().split('+').map(s => s.trim());
        const products = parts[2].trim().split('+').map(s => s.trim());

        // Parse each compound
        const allCompounds = [...reactants, ...products];
        const parsedCompounds = [];
        const allElements = new Set();

        for (let i = 0; i < allCompounds.length; i++) {
            const compound = allCompounds[i];
            if (!compound) {
                throw new Error('Empty compound found in equation');
            }

            // Remove any existing coefficients
            const cleanCompound = compound.replace(/^\d+/, '');
            if (!cleanCompound) {
                throw new Error(`Invalid compound: ${compound}`);
            }

            const elements = this.parseFormula(cleanCompound);
            if (Object.keys(elements).length === 0) {
                throw new Error(`Could not parse compound: ${compound}`);
            }

            parsedCompounds.push({
                formula: cleanCompound,
                elements: elements,
                isReactant: i < reactants.length
            });

            Object.keys(elements).forEach(el => allElements.add(el));
        }

        this.compounds = parsedCompounds;
        this.elements = Array.from(allElements).sort();

        return {
            reactants: reactants.map(r => r.replace(/^\d+/, '')),
            products: products.map(p => p.replace(/^\d+/, '')),
            compounds: parsedCompounds,
            elements: this.elements
        };
    }

    // Build the matrix for the system of equations
    buildMatrix() {
        const numCompounds = this.compounds.length;
        const numElements = this.elements.length;

        this.matrix = [];

        for (let i = 0; i < numElements; i++) {
            const row = [];
            const element = this.elements[i];

            for (let j = 0; j < numCompounds; j++) {
                const compound = this.compounds[j];
                const count = compound.elements[element] || 0;
                // Reactants are positive, products are negative
                row.push(compound.isReactant ? count : -count);
            }

            this.matrix.push(row);
        }

        return this.matrix;
    }

    // Solve using proper null space computation (following the pseudocode)
    solveMatrix() {
        const matrix = this.matrix.map(row => [...row]);
        const numRows = matrix.length;
        const numCols = matrix[0].length;

        // Helper functions
        const gcd = (a, b) => b === 0 ? Math.abs(a) : gcd(b, a % b);
        const lcm = (a, b) => Math.abs(a * b) / gcd(a, b);
        const Fraction = (num, den = 1) => {
            const g = gcd(Math.abs(num), Math.abs(den));
            return { num: num / g, den: den / g };
        };

        // Convert matrix to fractions for exact arithmetic
        const fracMatrix = matrix.map(row =>
            row.map(val => Fraction(Math.round(val * 1000), 1000))
        );

        // Perform RREF (Reduced Row Echelon Form)
        const rref = this.computeRREF(fracMatrix);

        // Compute null space from RREF
        const nullBasis = this.computeNullspaceFromRREF(rref);

        if (nullBasis.length === 0) {
            throw new Error('No nontrivial solution exists');
        }

        // Pick first basis vector (it's in fractions)
        const basisVector = nullBasis[0];

        // Find LCM of all denominators
        let commonDen = 1;
        basisVector.forEach(frac => {
            if (frac.den !== 1) {
                commonDen = lcm(commonDen, Math.abs(frac.den));
            }
        });

        // Multiply by LCM to get integers
        const integers = basisVector.map(frac => (frac.num * commonDen) / frac.den);

        // Ensure all positive (flip sign if needed)
        const hasNegative = integers.some(x => x < 0);
        if (hasNegative) {
            integers.forEach((val, i) => {
                integers[i] = -val;
            });
        }

        // Divide by GCD to get minimal integers
        const allGcd = integers.reduce((acc, val) =>
            val === 0 ? acc : gcd(acc, Math.abs(val)), Math.abs(integers.find(x => x !== 0) || 1)
        );

        const minimalInts = integers.map(val => Math.round(val / allGcd));

        // Ensure all coefficients are positive integers
        const finalSolution = minimalInts.map(x => Math.max(1, Math.abs(x)));

        this.solution = finalSolution;
        return finalSolution;
    }

    // Compute RREF using exact fraction arithmetic
    computeRREF(fracMatrix) {
        const matrix = fracMatrix.map(row => row.map(frac => ({ ...frac })));
        const numRows = matrix.length;
        const numCols = matrix[0].length;

        const addFractions = (a, b) => {
            const num = a.num * b.den + b.num * a.den;
            const den = a.den * b.den;
            const g = this.gcd(Math.abs(num), Math.abs(den));
            return { num: num / g, den: den / g };
        };

        const multiplyFractions = (a, b) => {
            const num = a.num * b.num;
            const den = a.den * b.den;
            const g = this.gcd(Math.abs(num), Math.abs(den));
            return { num: num / g, den: den / g };
        };

        const divideFractions = (a, b) => {
            if (b.num === 0) return { num: 0, den: 1 };
            return multiplyFractions(a, { num: b.den, den: b.num });
        };

        let currentRow = 0;

        for (let col = 0; col < numCols && currentRow < numRows; col++) {
            // Find pivot
            let pivotRow = -1;
            for (let row = currentRow; row < numRows; row++) {
                if (matrix[row][col].num !== 0) {
                    pivotRow = row;
                    break;
                }
            }

            if (pivotRow === -1) continue; // No pivot in this column

            // Swap rows
            if (pivotRow !== currentRow) {
                [matrix[currentRow], matrix[pivotRow]] = [matrix[pivotRow], matrix[currentRow]];
            }

            // Scale pivot row to make pivot = 1
            const pivot = matrix[currentRow][col];
            if (pivot.num !== 0) {
                for (let j = 0; j < numCols; j++) {
                    matrix[currentRow][j] = divideFractions(matrix[currentRow][j], pivot);
                }
            }

            // Eliminate column
            for (let row = 0; row < numRows; row++) {
                if (row !== currentRow && matrix[row][col].num !== 0) {
                    const factor = matrix[row][col];
                    for (let j = 0; j < numCols; j++) {
                        const term = multiplyFractions(factor, matrix[currentRow][j]);
                        matrix[row][j] = addFractions(matrix[row][j], { num: -term.num, den: term.den });
                    }
                }
            }

            currentRow++;
        }

        return matrix;
    }

    // Compute null space from RREF matrix
    computeNullspaceFromRREF(rref) {
        const numRows = rref.length;
        const numCols = rref[0].length;

        // Find pivot columns
        const pivotCols = [];
        const freeVars = [];

        for (let row = 0; row < numRows; row++) {
            for (let col = 0; col < numCols; col++) {
                if (rref[row][col].num !== 0) {
                    pivotCols.push(col);
                    break;
                }
            }
        }

        // Identify free variables
        for (let col = 0; col < numCols; col++) {
            if (!pivotCols.includes(col)) {
                freeVars.push(col);
            }
        }

        if (freeVars.length === 0) {
            return []; // No free variables, no null space
        }

        // Generate basis vectors for null space
        const nullBasis = [];

        for (const freeVar of freeVars) {
            const basisVector = new Array(numCols).fill(null).map(() => ({ num: 0, den: 1 }));
            basisVector[freeVar] = { num: 1, den: 1 }; // Set free variable to 1

            // Back substitute to find values of pivot variables
            for (let row = pivotCols.length - 1; row >= 0; row--) {
                const pivotCol = pivotCols[row];
                if (pivotCol < numCols) {
                    let sum = { num: 0, den: 1 };

                    for (let col = pivotCol + 1; col < numCols; col++) {
                        if (rref[row][col].num !== 0) {
                            const term = this.multiplyFractions(rref[row][col], basisVector[col]);
                            sum = this.addFractions(sum, term);
                        }
                    }

                    basisVector[pivotCol] = { num: -sum.num, den: sum.den };
                }
            }

            nullBasis.push(basisVector);
        }

        return nullBasis;
    }

    // Helper methods for fraction arithmetic
    gcd(a, b) {
        return b === 0 ? Math.abs(a) : this.gcd(b, a % b);
    }

    addFractions(a, b) {
        const num = a.num * b.den + b.num * a.den;
        const den = a.den * b.den;
        const g = this.gcd(Math.abs(num), Math.abs(den));
        return { num: num / g, den: den / g };
    }

    multiplyFractions(a, b) {
        const num = a.num * b.num;
        const den = a.den * b.den;
        const g = this.gcd(Math.abs(num), Math.abs(den));
        return { num: num / g, den: den / g };
    }

    // Balance the equation
    balance(equation) {
        try {
            const parsed = this.parseEquation(equation);
            this.buildMatrix();
            const coefficients = this.solveMatrix();

            // Verify the solution
            if (!this.verifySolution(coefficients)) {
                throw new Error('Could not find a valid solution');
            }

            return {
                success: true,
                coefficients: coefficients,
                compounds: this.compounds,
                elements: this.elements,
                matrix: this.matrix,
                balancedEquation: this.formatBalancedEquation(coefficients)
            };
        } catch (error) {
            return {
                success: false,
                error: error.message
            };
        }
    }

    // Verify that the solution balances the equation
    verifySolution(coefficients) {
        for (let i = 0; i < this.elements.length; i++) {
            let reactantSum = 0;
            let productSum = 0;

            for (let j = 0; j < this.compounds.length; j++) {
                const element = this.elements[i];
                const count = this.compounds[j].elements[element] || 0;
                const coeff = coefficients[j];

                if (this.compounds[j].isReactant) {
                    reactantSum += count * coeff;
                } else {
                    productSum += count * coeff;
                }
            }

            if (Math.abs(reactantSum - productSum) > 1e-10) {
                console.log(`Element ${this.elements[i]}: ${reactantSum} ≠ ${productSum}`);
                return false;
            }
        }
        return true;
    }

    // Format the balanced equation
    formatBalancedEquation(coefficients) {
        const reactants = [];
        const products = [];

        this.compounds.forEach((compound, i) => {
            const coeff = coefficients[i];
            const formula = compound.formula;
            const formatted = coeff === 1 ? formula : `${coeff}${formula}`;

            if (compound.isReactant) {
                reactants.push(formatted);
            } else {
                products.push(formatted);
            }
        });

        return `${reactants.join(' + ')} → ${products.join(' + ')}`;
    }
}

// Compound Database
const COMPOUND_DATABASE = [
    // Common gases
    { formula: "H₂", name: "Hydrogen gas", category: "Gas" },
    { formula: "O₂", name: "Oxygen gas", category: "Gas" },
    { formula: "N₂", name: "Nitrogen gas", category: "Gas" },
    { formula: "CO₂", name: "Carbon dioxide", category: "Gas" },
    { formula: "CO", name: "Carbon monoxide", category: "Gas" },
    { formula: "NH₃", name: "Ammonia", category: "Gas" },
    { formula: "CH₄", name: "Methane", category: "Gas" },
    { formula: "C₂H₆", name: "Ethane", category: "Gas" },
    { formula: "C₃H₈", name: "Propane", category: "Gas" },
    { formula: "C₄H₁₀", name: "Butane", category: "Gas" },
    { formula: "C₂H₄", name: "Ethylene", category: "Gas" },
    { formula: "C₂H₂", name: "Acetylene", category: "Gas" },
    { formula: "SO₂", name: "Sulfur dioxide", category: "Gas" },
    { formula: "SO₃", name: "Sulfur trioxide", category: "Gas" },
    { formula: "H₂S", name: "Hydrogen sulfide", category: "Gas" },
    { formula: "NO", name: "Nitric oxide", category: "Gas" },
    { formula: "NO₂", name: "Nitrogen dioxide", category: "Gas" },
    { formula: "N₂O", name: "Nitrous oxide", category: "Gas" },
    { formula: "N₂O₄", name: "Dinitrogen tetroxide", category: "Gas" },
    { formula: "Cl₂", name: "Chlorine gas", category: "Gas" },
    { formula: "F₂", name: "Fluorine gas", category: "Gas" },
    { formula: "Br₂", name: "Bromine gas", category: "Gas" },
    { formula: "I₂", name: "Iodine gas", category: "Gas" },
    { formula: "HCN", name: "Hydrogen cyanide", category: "Gas" },
    { formula: "PH₃", name: "Phosphine", category: "Gas" },

    // Water and common liquids
    { formula: "H₂O", name: "Water", category: "Liquid" },
    { formula: "H₂O₂", name: "Hydrogen peroxide", category: "Liquid" },
    { formula: "C₂H₅OH", name: "Ethanol", category: "Liquid" },
    { formula: "CH₃OH", name: "Methanol", category: "Liquid" },
    { formula: "C₃H₇OH", name: "Propanol", category: "Liquid" },
    { formula: "C₆H₆", name: "Benzene", category: "Liquid" },
    { formula: "CCl₄", name: "Carbon tetrachloride", category: "Liquid" },
    { formula: "CHCl₃", name: "Chloroform", category: "Liquid" },
    { formula: "CS₂", name: "Carbon disulfide", category: "Liquid" },

    // Common acids
    { formula: "HCl", name: "Hydrochloric acid", category: "Acid" },
    { formula: "H₂SO₄", name: "Sulfuric acid", category: "Acid" },
    { formula: "HNO₃", name: "Nitric acid", category: "Acid" },
    { formula: "CH₃COOH", name: "Acetic acid", category: "Acid" },
    { formula: "H₃PO₄", name: "Phosphoric acid", category: "Acid" },
    { formula: "H₂CO₃", name: "Carbonic acid", category: "Acid" },
    { formula: "HF", name: "Hydrofluoric acid", category: "Acid" },
    { formula: "HBr", name: "Hydrobromic acid", category: "Acid" },
    { formula: "HI", name: "Hydroiodic acid", category: "Acid" },
    { formula: "H₂S", name: "Hydrosulfuric acid", category: "Acid" },
    { formula: "HClO₄", name: "Perchloric acid", category: "Acid" },
    { formula: "HClO₃", name: "Chloric acid", category: "Acid" },
    { formula: "HClO₂", name: "Chlorous acid", category: "Acid" },
    { formula: "HClO", name: "Hypochlorous acid", category: "Acid" },

    // Common bases
    { formula: "NaOH", name: "Sodium hydroxide", category: "Base" },
    { formula: "KOH", name: "Potassium hydroxide", category: "Base" },
    { formula: "Ca(OH)₂", name: "Calcium hydroxide", category: "Base" },
    { formula: "Mg(OH)₂", name: "Magnesium hydroxide", category: "Base" },
    { formula: "Ba(OH)₂", name: "Barium hydroxide", category: "Base" },
    { formula: "Al(OH)₃", name: "Aluminum hydroxide", category: "Base" },
    { formula: "Fe(OH)₂", name: "Iron(II) hydroxide", category: "Base" },
    { formula: "Fe(OH)₃", name: "Iron(III) hydroxide", category: "Base" },
    { formula: "Cu(OH)₂", name: "Copper(II) hydroxide", category: "Base" },
    { formula: "Zn(OH)₂", name: "Zinc hydroxide", category: "Base" },
    { formula: "LiOH", name: "Lithium hydroxide", category: "Base" },
    { formula: "RbOH", name: "Rubidium hydroxide", category: "Base" },
    { formula: "CsOH", name: "Cesium hydroxide", category: "Base" },

    // Sulfates (including Iron(III) sulfate!)
    { formula: "Fe₂(SO₄)₃", name: "Iron(III) sulfate", category: "Salt" },
    { formula: "FeSO₄", name: "Iron(II) sulfate", category: "Salt" },
    { formula: "CuSO₄", name: "Copper(II) sulfate", category: "Salt" },
    { formula: "ZnSO₄", name: "Zinc sulfate", category: "Salt" },
    { formula: "Al₂(SO₄)₃", name: "Aluminum sulfate", category: "Salt" },
    { formula: "Na₂SO₄", name: "Sodium sulfate", category: "Salt" },
    { formula: "K₂SO₄", name: "Potassium sulfate", category: "Salt" },
    { formula: "CaSO₄", name: "Calcium sulfate", category: "Salt" },
    { formula: "MgSO₄", name: "Magnesium sulfate", category: "Salt" },
    { formula: "BaSO₄", name: "Barium sulfate", category: "Salt" },
    { formula: "PbSO₄", name: "Lead(II) sulfate", category: "Salt" },
    { formula: "Pb(SO₄)₂", name: "Lead(IV) sulfate", category: "Salt" },
    { formula: "Ag₂SO₄", name: "Silver sulfate", category: "Salt" },

    // Nitrates
    { formula: "AgNO₃", name: "Silver nitrate", category: "Salt" },
    { formula: "Pb(NO₃)₂", name: "Lead(II) nitrate", category: "Salt" },
    { formula: "Pb(NO₃)₄", name: "Lead(IV) nitrate", category: "Salt" },
    { formula: "Cu(NO₃)₂", name: "Copper(II) nitrate", category: "Salt" },
    { formula: "Fe(NO₃)₃", name: "Iron(III) nitrate", category: "Salt" },
    { formula: "Fe(NO₃)₂", name: "Iron(II) nitrate", category: "Salt" },
    { formula: "Al(NO₃)₃", name: "Aluminum nitrate", category: "Salt" },
    { formula: "Zn(NO₃)₂", name: "Zinc nitrate", category: "Salt" },
    { formula: "Ca(NO₃)₂", name: "Calcium nitrate", category: "Salt" },
    { formula: "Mg(NO₃)₂", name: "Magnesium nitrate", category: "Salt" },
    { formula: "Ba(NO₃)₂", name: "Barium nitrate", category: "Salt" },
    { formula: "NaNO₃", name: "Sodium nitrate", category: "Salt" },
    { formula: "KNO₃", name: "Potassium nitrate", category: "Salt" },
    { formula: "NH₄NO₃", name: "Ammonium nitrate", category: "Salt" },

    // Chlorides
    { formula: "NaCl", name: "Sodium chloride", category: "Salt" },
    { formula: "KCl", name: "Potassium chloride", category: "Salt" },
    { formula: "CaCl₂", name: "Calcium chloride", category: "Salt" },
    { formula: "MgCl₂", name: "Magnesium chloride", category: "Salt" },
    { formula: "FeCl₃", name: "Iron(III) chloride", category: "Salt" },
    { formula: "FeCl₂", name: "Iron(II) chloride", category: "Salt" },
    { formula: "CuCl₂", name: "Copper(II) chloride", category: "Salt" },
    { formula: "CuCl", name: "Copper(I) chloride", category: "Salt" },
    { formula: "ZnCl₂", name: "Zinc chloride", category: "Salt" },
    { formula: "AlCl₃", name: "Aluminum chloride", category: "Salt" },
    { formula: "AgCl", name: "Silver chloride", category: "Salt" },
    { formula: "PbCl₂", name: "Lead(II) chloride", category: "Salt" },
    { formula: "PbCl₄", name: "Lead(IV) chloride", category: "Salt" },
    { formula: "BaCl₂", name: "Barium chloride", category: "Salt" },
    { formula: "NH₄Cl", name: "Ammonium chloride", category: "Salt" },

    // Bromides
    { formula: "NaBr", name: "Sodium bromide", category: "Salt" },
    { formula: "KBr", name: "Potassium bromide", category: "Salt" },
    { formula: "CaBr₂", name: "Calcium bromide", category: "Salt" },
    { formula: "MgBr₂", name: "Magnesium bromide", category: "Salt" },
    { formula: "AgBr", name: "Silver bromide", category: "Salt" },

    // Iodides
    { formula: "NaI", name: "Sodium iodide", category: "Salt" },
    { formula: "KI", name: "Potassium iodide", category: "Salt" },
    { formula: "CaI₂", name: "Calcium iodide", category: "Salt" },
    { formula: "MgI₂", name: "Magnesium iodide", category: "Salt" },
    { formula: "AgI", name: "Silver iodide", category: "Salt" },
    { formula: "PbI₂", name: "Lead(II) iodide", category: "Salt" },

    // Carbonates and Bicarbonates
    { formula: "Na₂CO₃", name: "Sodium carbonate", category: "Salt" },
    { formula: "K₂CO₃", name: "Potassium carbonate", category: "Salt" },
    { formula: "CaCO₃", name: "Calcium carbonate", category: "Salt" },
    { formula: "MgCO₃", name: "Magnesium carbonate", category: "Salt" },
    { formula: "BaCO₃", name: "Barium carbonate", category: "Salt" },
    { formula: "FeCO₃", name: "Iron(II) carbonate", category: "Salt" },
    { formula: "CuCO₃", name: "Copper(II) carbonate", category: "Salt" },
    { formula: "ZnCO₃", name: "Zinc carbonate", category: "Salt" },
    { formula: "NaHCO₃", name: "Sodium bicarbonate", category: "Salt" },
    { formula: "KHCO₃", name: "Potassium bicarbonate", category: "Salt" },
    { formula: "Ca(HCO₃)₂", name: "Calcium bicarbonate", category: "Salt" },
    { formula: "Mg(HCO₃)₂", name: "Magnesium bicarbonate", category: "Salt" },

    // Phosphates
    { formula: "Na₃PO₄", name: "Sodium phosphate", category: "Salt" },
    { formula: "K₃PO₄", name: "Potassium phosphate", category: "Salt" },
    { formula: "Ca₃(PO₄)₂", name: "Calcium phosphate", category: "Salt" },
    { formula: "Mg₃(PO₄)₂", name: "Magnesium phosphate", category: "Salt" },
    { formula: "FePO₄", name: "Iron(III) phosphate", category: "Salt" },
    { formula: "Fe₃(PO₄)₂", name: "Iron(II) phosphate", category: "Salt" },
    { formula: "AlPO₄", name: "Aluminum phosphate", category: "Salt" },
    { formula: "Pb₃(PO₄)₂", name: "Lead(II) phosphate", category: "Salt" },
    { formula: "Pb₃(PO₄)₄", name: "Lead(IV) phosphate", category: "Salt" },
    { formula: "Cu₃(PO₄)₂", name: "Copper(II) phosphate", category: "Salt" },
    { formula: "Zn₃(PO₄)₂", name: "Zinc phosphate", category: "Salt" },
    { formula: "Ag₃PO₄", name: "Silver phosphate", category: "Salt" },
    { formula: "Ba₃(PO₄)₂", name: "Barium phosphate", category: "Salt" },
    { formula: "Sr₃(PO₄)₂", name: "Strontium phosphate", category: "Salt" },
    { formula: "Mn₃(PO₄)₂", name: "Manganese(II) phosphate", category: "Salt" },
    { formula: "Co₃(PO₄)₂", name: "Cobalt(II) phosphate", category: "Salt" },
    { formula: "Ni₃(PO₄)₂", name: "Nickel(II) phosphate", category: "Salt" },
    { formula: "Cr₃(PO₄)₄", name: "Chromium(III) phosphate", category: "Salt" },
    { formula: "(NH₄)₃PO₄", name: "Ammonium phosphate", category: "Salt" },
    { formula: "(NH₄)₂HPO₄", name: "Diammonium phosphate", category: "Salt" },
    { formula: "NH₄H₂PO₄", name: "Monoammonium phosphate", category: "Salt" },

    // Ammonium compounds
    { formula: "NH₄Cl", name: "Ammonium chloride", category: "Salt" },
    { formula: "(NH₄)₂SO₄", name: "Ammonium sulfate", category: "Salt" },
    { formula: "NH₄NO₃", name: "Ammonium nitrate", category: "Salt" },
    { formula: "(NH₄)₃PO₄", name: "Ammonium phosphate", category: "Salt" },
    { formula: "(NH₄)₂CO₃", name: "Ammonium carbonate", category: "Salt" },
    { formula: "NH₄OH", name: "Ammonium hydroxide", category: "Base" },

    // Common oxides
    { formula: "Fe₂O₃", name: "Iron(III) oxide", category: "Oxide" },
    { formula: "FeO", name: "Iron(II) oxide", category: "Oxide" },
    { formula: "Fe₃O₄", name: "Iron(II,III) oxide", category: "Oxide" },
    { formula: "Al₂O₃", name: "Aluminum oxide", category: "Oxide" },
    { formula: "CuO", name: "Copper(II) oxide", category: "Oxide" },
    { formula: "Cu₂O", name: "Copper(I) oxide", category: "Oxide" },
    { formula: "ZnO", name: "Zinc oxide", category: "Oxide" },
    { formula: "MgO", name: "Magnesium oxide", category: "Oxide" },
    { formula: "CaO", name: "Calcium oxide", category: "Oxide" },
    { formula: "BaO", name: "Barium oxide", category: "Oxide" },
    { formula: "Na₂O", name: "Sodium oxide", category: "Oxide" },
    { formula: "K₂O", name: "Potassium oxide", category: "Oxide" },
    { formula: "Li₂O", name: "Lithium oxide", category: "Oxide" },
    { formula: "PbO", name: "Lead(II) oxide", category: "Oxide" },
    { formula: "PbO₂", name: "Lead(IV) oxide", category: "Oxide" },
    { formula: "Ag₂O", name: "Silver oxide", category: "Oxide" },
    { formula: "MnO", name: "Manganese(II) oxide", category: "Oxide" },
    { formula: "MnO₂", name: "Manganese(IV) oxide", category: "Oxide" },
    { formula: "Mn₂O₃", name: "Manganese(III) oxide", category: "Oxide" },
    { formula: "Cr₂O₃", name: "Chromium(III) oxide", category: "Oxide" },
    { formula: "CrO₃", name: "Chromium(VI) oxide", category: "Oxide" },

    // Elements (metals)
    { formula: "Fe", name: "Iron", category: "Element" },
    { formula: "Cu", name: "Copper", category: "Element" },
    { formula: "Zn", name: "Zinc", category: "Element" },
    { formula: "Al", name: "Aluminum", category: "Element" },
    { formula: "Mg", name: "Magnesium", category: "Element" },
    { formula: "Ca", name: "Calcium", category: "Element" },
    { formula: "Na", name: "Sodium", category: "Element" },
    { formula: "K", name: "Potassium", category: "Element" },
    { formula: "Li", name: "Lithium", category: "Element" },
    { formula: "Ag", name: "Silver", category: "Element" },
    { formula: "Pb", name: "Lead", category: "Element" },
    { formula: "Sn", name: "Tin", category: "Element" },
    { formula: "Ni", name: "Nickel", category: "Element" },
    { formula: "Co", name: "Cobalt", category: "Element" },
    { formula: "Mn", name: "Manganese", category: "Element" },
    { formula: "Cr", name: "Chromium", category: "Element" },
    { formula: "Ba", name: "Barium", category: "Element" },
    { formula: "Sr", name: "Strontium", category: "Element" },

    // Elements (nonmetals)
    { formula: "C", name: "Carbon", category: "Element" },
    { formula: "S", name: "Sulfur", category: "Element" },
    { formula: "P", name: "Phosphorus", category: "Element" },
    { formula: "Si", name: "Silicon", category: "Element" },
    { formula: "B", name: "Boron", category: "Element" },

    // Organic compounds
    { formula: "C₆H₁₂O₆", name: "Glucose", category: "Organic" },
    { formula: "C₁₂H₂₂O₁₁", name: "Sucrose", category: "Organic" },
    { formula: "C₂H₄O₂", name: "Acetic acid", category: "Organic" },
    { formula: "C₃H₆O₃", name: "Lactic acid", category: "Organic" },
    { formula: "C₄H₆O₄", name: "Succinic acid", category: "Organic" },
    { formula: "C₆H₈O₇", name: "Citric acid", category: "Organic" },
    { formula: "C₈H₁₈", name: "Octane", category: "Organic" },
    { formula: "C₁₀H₈", name: "Naphthalene", category: "Organic" },

    // Peroxides and special compounds
    { formula: "Na₂O₂", name: "Sodium peroxide", category: "Peroxide" },
    { formula: "BaO₂", name: "Barium peroxide", category: "Peroxide" },
    { formula: "CaC₂", name: "Calcium carbide", category: "Carbide" },
    { formula: "Al₄C₃", name: "Aluminum carbide", category: "Carbide" },
    { formula: "SiC", name: "Silicon carbide", category: "Carbide" },
    { formula: "CaF₂", name: "Calcium fluoride", category: "Fluoride" },
    { formula: "NaF", name: "Sodium fluoride", category: "Fluoride" },
    { formula: "AlF₃", name: "Aluminum fluoride", category: "Fluoride" },

    // Additional Lead Compounds
    { formula: "PbO", name: "Lead(II) oxide", category: "Oxide" },
    { formula: "PbO₂", name: "Lead(IV) oxide", category: "Oxide" },
    { formula: "Pb₃O₄", name: "Lead(II,IV) oxide", category: "Oxide" },
    { formula: "PbS", name: "Lead(II) sulfide", category: "Sulfide" },
    { formula: "PbCO₃", name: "Lead(II) carbonate", category: "Salt" },
    { formula: "Pb(OH)₂", name: "Lead(II) hydroxide", category: "Base" },
    { formula: "Pb(OH)₄", name: "Lead(IV) hydroxide", category: "Base" },

    // More Tin Compounds
    { formula: "SnO", name: "Tin(II) oxide", category: "Oxide" },
    { formula: "SnO₂", name: "Tin(IV) oxide", category: "Oxide" },
    { formula: "SnCl₂", name: "Tin(II) chloride", category: "Salt" },
    { formula: "SnCl₄", name: "Tin(IV) chloride", category: "Salt" },
    { formula: "Sn(SO₄)₂", name: "Tin(IV) sulfate", category: "Salt" },
    { formula: "SnSO₄", name: "Tin(II) sulfate", category: "Salt" },

    // More Chromium Compounds
    { formula: "CrCl₂", name: "Chromium(II) chloride", category: "Salt" },
    { formula: "CrCl₃", name: "Chromium(III) chloride", category: "Salt" },
    { formula: "Cr₂(SO₄)₃", name: "Chromium(III) sulfate", category: "Salt" },
    { formula: "CrSO₄", name: "Chromium(II) sulfate", category: "Salt" },
    { formula: "K₂Cr₂O₇", name: "Potassium dichromate", category: "Salt" },
    { formula: "K₂CrO₄", name: "Potassium chromate", category: "Salt" },

    // More Manganese Compounds
    { formula: "MnCl₂", name: "Manganese(II) chloride", category: "Salt" },
    { formula: "MnSO₄", name: "Manganese(II) sulfate", category: "Salt" },
    { formula: "Mn(NO₃)₂", name: "Manganese(II) nitrate", category: "Salt" },
    { formula: "KMnO₄", name: "Potassium permanganate", category: "Salt" },

    // More Cobalt and Nickel Compounds
    { formula: "CoCl₂", name: "Cobalt(II) chloride", category: "Salt" },
    { formula: "CoSO₄", name: "Cobalt(II) sulfate", category: "Salt" },
    { formula: "Co(NO₃)₂", name: "Cobalt(II) nitrate", category: "Salt" },
    { formula: "NiCl₂", name: "Nickel(II) chloride", category: "Salt" },
    { formula: "NiSO₄", name: "Nickel(II) sulfate", category: "Salt" },
    { formula: "Ni(NO₃)₂", name: "Nickel(II) nitrate", category: "Salt" },

    // Sulfides
    { formula: "H₂S", name: "Hydrogen sulfide", category: "Sulfide" },
    { formula: "Na₂S", name: "Sodium sulfide", category: "Sulfide" },
    { formula: "K₂S", name: "Potassium sulfide", category: "Sulfide" },
    { formula: "CaS", name: "Calcium sulfide", category: "Sulfide" },
    { formula: "MgS", name: "Magnesium sulfide", category: "Sulfide" },
    { formula: "FeS", name: "Iron(II) sulfide", category: "Sulfide" },
    { formula: "Fe₂S₃", name: "Iron(III) sulfide", category: "Sulfide" },
    { formula: "CuS", name: "Copper(II) sulfide", category: "Sulfide" },
    { formula: "Cu₂S", name: "Copper(I) sulfide", category: "Sulfide" },
    { formula: "ZnS", name: "Zinc sulfide", category: "Sulfide" },
    { formula: "Ag₂S", name: "Silver sulfide", category: "Sulfide" },
    { formula: "Al₂S₃", name: "Aluminum sulfide", category: "Sulfide" },

    // More Complex Compounds
    { formula: "Ca(ClO)₂", name: "Calcium hypochlorite", category: "Salt" },
    { formula: "NaClO", name: "Sodium hypochlorite", category: "Salt" },
    { formula: "NaClO₂", name: "Sodium chlorite", category: "Salt" },
    { formula: "NaClO₃", name: "Sodium chlorate", category: "Salt" },
    { formula: "NaClO₄", name: "Sodium perchlorate", category: "Salt" },
    { formula: "KClO₃", name: "Potassium chlorate", category: "Salt" },
    { formula: "KClO₄", name: "Potassium perchlorate", category: "Salt" },

    // Hydrates (common in chemistry)
    { formula: "CuSO₄·5H₂O", name: "Copper(II) sulfate pentahydrate", category: "Hydrate" },
    { formula: "FeSO₄·7H₂O", name: "Iron(II) sulfate heptahydrate", category: "Hydrate" },
    { formula: "MgSO₄·7H₂O", name: "Magnesium sulfate heptahydrate", category: "Hydrate" },
    { formula: "CaCl₂·2H₂O", name: "Calcium chloride dihydrate", category: "Hydrate" },
    { formula: "Na₂CO₃·10H₂O", name: "Sodium carbonate decahydrate", category: "Hydrate" }
];

// UI Controller
class EquationBalancerUI {
    constructor() {
        this.balancer = new ChemicalEquationBalancer();
        this.compounds = COMPOUND_DATABASE;
        this.activeDropdown = null;
        this.showingSmartSuggestions = false;
        this.initializeEventListeners();
        this.setupDropdowns();
    }

    initializeEventListeners() {
        console.log('Initializing event listeners...');

        // Get elements
        const balanceBtn = document.getElementById('balance-btn');
        const clearBtn = document.getElementById('clear-btn');
        const reactantsInput = document.getElementById('reactants-input');
        const productsInput = document.getElementById('products-input');
        const copyBtn = document.getElementById('copy-btn');

        // Balance button
        if (balanceBtn) {
            console.log('Balance button found, adding event listener');
            balanceBtn.addEventListener('click', (e) => {
                e.preventDefault();
                console.log('Balance button clicked');
                this.balanceEquation();
            });
        } else {
            console.error('Balance button not found!');
        }

        // Clear button
        if (clearBtn) {
            clearBtn.addEventListener('click', (e) => {
                e.preventDefault();
                this.clearInputs();
            });
        }

        // Enter key support
        if (reactantsInput && productsInput) {
            [reactantsInput, productsInput].forEach(input => {
                input.addEventListener('keypress', (e) => {
                    if (e.key === 'Enter') {
                        e.preventDefault();
                        this.balanceEquation();
                    }
                });
            });
        }

        // Copy button
        if (copyBtn) {
            copyBtn.addEventListener('click', () => this.copyResult());
        }

        // Add example placeholders
        this.addExamplePlaceholders();

        // Setup dropdown functionality
        this.setupDropdownEvents();

        // Close dropdowns when clicking outside
        document.addEventListener('click', (e) => {
            if (!e.target.closest('.relative')) {
                this.closeAllDropdowns();
            }
        });
    }

    setupDropdowns() {
        this.populateDropdown('reactants');
        this.populateDropdown('products');
    }

    setupDropdownEvents() {
        // Reactants dropdown
        const reactantsBtn = document.getElementById('reactants-dropdown-btn');
        const reactantsDropdown = document.getElementById('reactants-dropdown');
        const reactantsInput = document.getElementById('reactants-input');

        if (reactantsBtn && reactantsDropdown) {
            reactantsBtn.addEventListener('click', (e) => {
                e.stopPropagation();
                this.toggleDropdown('reactants');
            });

            reactantsInput.addEventListener('focus', () => {
                console.log('Reactants input focused, value:', reactantsInput.value);

                // If the input ends with " + ", show smart suggestions
                if (reactantsInput.value.endsWith(' + ')) {
                    const compounds = reactantsInput.value.replace(' + ', '').split('+').map(c => c.trim()).filter(c => c);
                    if (compounds.length > 0) {
                        console.log('Showing smart suggestions on focus for:', compounds[compounds.length - 1]);
                        this.showSmartSuggestions('reactants', compounds[compounds.length - 1]);
                        return;
                    }
                }

                // Check if products field has content - show smart reactant suggestions
                const productsValue = productsInput.value.trim();
                if (productsValue && !reactantsInput.value.trim()) {
                    console.log('Products field has content, showing smart reactant suggestions');
                    this.showSmartSuggestionsBasedOnOtherField('reactants', productsValue);
                    return;
                }

                this.showDropdown('reactants');
                // Filter based on current input
                const query = this.getCurrentCompound(reactantsInput.value);
                this.filterDropdown('reactants', query);
            });

            // Handle space bar for smart suggestions
            reactantsInput.addEventListener('keydown', (e) => {
                if (e.key === ' ' && e.target.value.trim() && !e.target.value.endsWith(' ') && !e.target.value.endsWith('+')) {
                    console.log('Space bar pressed, current value:', e.target.value);
                    e.preventDefault(); // Prevent the space from being added normally
                    e.stopPropagation(); // Prevent other handlers from interfering

                    const currentValue = e.target.value.trim();
                    console.log('Adding + and showing suggestions for:', currentValue);

                    // Add " + " to the input
                    e.target.value = currentValue + ' + ';
                    e.target.setSelectionRange(e.target.value.length, e.target.value.length);

                    // Show smart suggestions with a small delay to ensure it's not immediately closed
                    setTimeout(() => {
                        console.log('Calling showSmartSuggestions...');
                        this.showSmartSuggestions('reactants', currentValue);

                        // Update product suggestions
                        this.updateProductSuggestions();
                    }, 10);
                }
            });

            reactantsInput.addEventListener('input', (e) => {
                console.log('Input event fired, value:', e.target.value, 'inputType:', e.inputType, 'showingSmartSuggestions:', this.showingSmartSuggestions);

                // Don't process if we're showing smart suggestions or value ends with " + "
                if (this.showingSmartSuggestions || e.target.value.endsWith(' + ')) {
                    console.log('Skipping input event - smart suggestions active or value ends with " + "');
                    return;
                }

                const query = this.getCurrentCompound(e.target.value);
                console.log('Processing input event normally, query:', query);

                this.showDropdown('reactants');
                this.filterDropdown('reactants', query);

                // Update product suggestions when reactants change
                this.updateProductSuggestions();

                // Update product placeholder in real-time
                this.updateProductPlaceholder();
            });
        }

        // Products dropdown
        const productsBtn = document.getElementById('products-dropdown-btn');
        const productsDropdown = document.getElementById('products-dropdown');
        const productsInput = document.getElementById('products-input');

        if (productsBtn && productsDropdown) {
            productsBtn.addEventListener('click', (e) => {
                e.stopPropagation();
                this.toggleDropdown('products');
            });

            productsInput.addEventListener('focus', () => {
                // Clear placeholder suggestions when focused
                if (productsInput.classList.contains('has-suggestion')) {
                    productsInput.value = '';
                    productsInput.classList.remove('has-suggestion');
                    productsInput.style.color = '';
                }

                // Check if reactants field has content - show smart product suggestions
                const reactantsValue = reactantsInput.value.trim();
                if (reactantsValue && !productsInput.value.trim()) {
                    console.log('Reactants field has content, showing smart product suggestions');
                    this.showSmartSuggestionsBasedOnOtherField('products', reactantsValue);
                    return;
                }

                this.showDropdown('products');
                const query = this.getCurrentCompound(productsInput.value);
                this.filterDropdown('products', query);
            });

            // Handle space bar and Tab key for products
            productsInput.addEventListener('keydown', (e) => {
                if (e.key === 'Tab' && productsInput.classList.contains('has-suggestion')) {
                    e.preventDefault();
                    productsInput.classList.remove('has-suggestion');
                    productsInput.style.color = '';
                } else if (e.key === ' ' && e.target.value.trim() && !e.target.value.endsWith(' ') && !e.target.value.endsWith('+')) {
                    e.preventDefault(); // Prevent the space from being added normally

                    const currentValue = e.target.value.trim();
                    // Add " + " to the input
                    e.target.value = currentValue + ' + ';
                    e.target.setSelectionRange(e.target.value.length, e.target.value.length);

                    // Show smart suggestions based on what was typed
                    this.showSmartSuggestions('products', currentValue);
                }
            });

            productsInput.addEventListener('input', (e) => {
                // Remove suggestion styling when user types
                if (productsInput.classList.contains('has-suggestion')) {
                    productsInput.classList.remove('has-suggestion');
                    productsInput.style.color = '';
                }

                const query = this.getCurrentCompound(e.target.value);

                // Don't process if we just handled a space bar event
                if (e.target.value.endsWith(' + ') && e.inputType !== 'insertText') {
                    return;
                }

                this.showDropdown('products');
                this.filterDropdown('products', query);

                // Update reactant placeholder in real-time
                this.updateReactantPlaceholder();
            });
        }
    }

    // Get the current compound being typed (after the last +)
    getCurrentCompound(inputValue) {
        return inputValue.split('+').pop().trim();
    }

    // Update product placeholder based on current reactants
    updateProductPlaceholder() {
        const reactantsInput = document.getElementById('reactants-input');
        const productsInput = document.getElementById('products-input');

        if (!reactantsInput || !productsInput) return;
        if (productsInput.value.trim()) return; // Don't override if user has typed something

        const reactants = reactantsInput.value.split('+').map(r => r.trim()).filter(r => r);
        console.log('updateProductPlaceholder - reactants:', reactants);

        if (reactants.length > 0) {
            const suggestedProducts = this.predictProducts(reactants);
            console.log('updateProductPlaceholder - suggestedProducts:', suggestedProducts);

            if (suggestedProducts.length > 0) {
                productsInput.placeholder = suggestedProducts.join(' + ');
                console.log('Set placeholder to:', suggestedProducts.join(' + '));
            } else {
                productsInput.placeholder = 'Products will appear here';
                console.log('No products found, set default placeholder');
            }
        } else {
            productsInput.placeholder = 'CO₂ + H₂O'; // Default placeholder
        }
    }

    // Update reactant placeholder based on current products
    updateReactantPlaceholder() {
        const reactantsInput = document.getElementById('reactants-input');
        const productsInput = document.getElementById('products-input');

        if (!reactantsInput || !productsInput) return;
        if (reactantsInput.value.trim()) return; // Don't override if user has typed something

        const products = productsInput.value.split('+').map(p => p.trim()).filter(p => p);
        if (products.length > 0) {
            const suggestedReactants = this.predictReactants(products);
            if (suggestedReactants.length > 0) {
                reactantsInput.placeholder = suggestedReactants.join(' + ');
            } else {
                reactantsInput.placeholder = 'Reactants will appear here';
            }
        } else {
            reactantsInput.placeholder = 'CH₄ + O₂'; // Default placeholder
        }
    }

    // Smart suggestions based on existing compounds
    getSmartSuggestions(type, currentCompounds) {
        console.log('getSmartSuggestions called with:', type, currentCompounds);
        const existing = currentCompounds.split('+').map(c => c.trim()).filter(c => c);
        console.log('Existing compounds:', existing);

        if (type === 'reactants') {
            const suggestions = this.getReactantSuggestions(existing);
            console.log('Reactant suggestions:', suggestions);
            return suggestions;
        } else {
            return this.getProductSuggestions(existing);
        }
    }

    // Get smart reactant suggestions based on what's already entered
    getReactantSuggestions(existingReactants) {
        const suggestions = [];

        // Common reaction patterns - much more comprehensive
        const patterns = {
            // Acids - suggest bases and metals
            'HCl': ['NaOH', 'KOH', 'Ca(OH)₂', 'Mg(OH)₂', 'Zn', 'Fe', 'Al', 'Mg', 'CaCO₃', 'Na₂CO₃'],
            'H₂SO₄': ['NaOH', 'KOH', 'Ca(OH)₂', 'Mg(OH)₂', 'Zn', 'Fe', 'Al', 'Mg', 'CaCO₃', 'Na₂CO₃'],
            'HNO₃': ['NaOH', 'KOH', 'Ca(OH)₂', 'Mg(OH)₂', 'Cu', 'Ag', 'Zn', 'Fe', 'Al'],
            'CH₃COOH': ['NaOH', 'KOH', 'Ca(OH)₂', 'Mg(OH)₂', 'Na₂CO₃', 'CaCO₃'],
            'H₃PO₄': ['NaOH', 'KOH', 'Ca(OH)₂', 'Mg(OH)₂'],

            // Bases - suggest acids
            'NaOH': ['HCl', 'H₂SO₄', 'HNO₃', 'CH₃COOH', 'H₃PO₄', 'CO₂'],
            'KOH': ['HCl', 'H₂SO₄', 'HNO₃', 'CH₃COOH', 'H₃PO₄', 'CO₂'],
            'Ca(OH)₂': ['HCl', 'H₂SO₄', 'HNO₃', 'CH₃COOH', 'CO₂'],
            'Mg(OH)₂': ['HCl', 'H₂SO₄', 'HNO₃', 'CH₃COOH'],

            // Metals - suggest oxygen, acids, water, or salts
            'Fe': ['O₂', 'HCl', 'H₂SO₄', 'CuSO₄', 'AgNO₃', 'H₂O', 'S'],
            'Cu': ['O₂', 'HNO₃', 'H₂SO₄', 'AgNO₃', 'S'],
            'Zn': ['O₂', 'HCl', 'H₂SO₄', 'CuSO₄', 'AgNO₃', 'H₂O'],
            'Al': ['O₂', 'HCl', 'H₂SO₄', 'CuSO₄', 'AgNO₃', 'Fe₂O₃'],
            'Mg': ['O₂', 'HCl', 'H₂SO₄', 'CuSO₄', 'AgNO₃', 'H₂O'],
            'Ca': ['O₂', 'HCl', 'H₂SO₄', 'H₂O'],
            'Na': ['O₂', 'H₂O', 'Cl₂'],
            'K': ['O₂', 'H₂O', 'Cl₂'],

            // Hydrocarbons - suggest oxygen for combustion
            'CH₄': ['O₂'],
            'C₂H₆': ['O₂'],
            'C₃H₈': ['O₂'],
            'C₄H₁₀': ['O₂'],
            'C₂H₄': ['O₂', 'Br₂', 'H₂O'],
            'C₂H₂': ['O₂', 'H₂O'],
            'C₆H₆': ['O₂', 'Br₂'],
            'C₂H₅OH': ['O₂', 'K₂Cr₂O₇', 'KMnO₄'],
            'CH₃OH': ['O₂', 'K₂Cr₂O₇', 'KMnO₄'],
            'C₆H₁₂O₆': ['O₂', 'C₂H₅OH'],

            // Oxygen - suggest fuels and metals
            'O₂': ['CH₄', 'C₂H₆', 'C₃H₈', 'C₂H₄', 'C₂H₂', 'Fe', 'Cu', 'Zn', 'Al', 'Mg', 'Ca', 'Na', 'K', 'S', 'P'],

            // Salts for double displacement
            'NaCl': ['AgNO₃', 'Pb(NO₃)₂', 'Ca(NO₃)₂', 'H₂SO₄'],
            'AgNO₃': ['NaCl', 'KCl', 'CaCl₂', 'HCl'],
            'Pb(NO₃)₂': ['NaCl', 'KCl', 'CaCl₂', 'H₂SO₄'],
            'CuSO₄': ['NaOH', 'KOH', 'BaCl₂', 'Zn', 'Fe', 'Al'],
            'FeCl₃': ['NaOH', 'KOH', 'AgNO₃'],

            // More salts
            'BaCl₂': ['H₂SO₄', 'Na₂SO₄', 'K₂SO₄', 'AgNO₃'],
            'MgSO₄': ['NaOH', 'KOH', 'BaCl₂'],
            'ZnSO₄': ['NaOH', 'KOH', 'BaCl₂', 'Fe', 'Al'],
            'FeSO₄': ['NaOH', 'KOH', 'BaCl₂'],
            'AlCl₃': ['NaOH', 'KOH', 'AgNO₃'],
            'CaCl₂': ['H₂SO₄', 'Na₂SO₄', 'AgNO₃', 'Na₂CO₃'],
            'MgCl₂': ['H₂SO₄', 'Na₂SO₄', 'AgNO₃', 'NaOH'],

            // Carbonates - suggest acids
            'CaCO₃': ['HCl', 'H₂SO₄', 'HNO₃', 'CH₃COOH'],
            'Na₂CO₃': ['HCl', 'H₂SO₄', 'HNO₃', 'CH₃COOH'],
            'MgCO₃': ['HCl', 'H₂SO₄', 'HNO₃'],

            // Oxides - suggest acids or water
            'CaO': ['H₂O', 'HCl', 'H₂SO₄', 'CO₂'],
            'MgO': ['H₂O', 'HCl', 'H₂SO₄'],
            'Fe₂O₃': ['HCl', 'H₂SO₄', 'Al', 'CO'],
            'CuO': ['HCl', 'H₂SO₄', 'H₂'],
            'Al₂O₃': ['HCl', 'H₂SO₄'],

            // Gases
            'H₂': ['O₂', 'Cl₂', 'Br₂', 'CuO', 'Fe₂O₃'],
            'Cl₂': ['H₂', 'Na', 'K', 'Fe', 'NaBr', 'KI'],
            'Br₂': ['H₂', 'KI', 'NaI', 'C₂H₄', 'C₆H₆'],
            'I₂': ['H₂', 'Na₂S₂O₃'],
            'F₂': ['H₂', 'NaCl', 'KCl', 'KBr', 'KI'],
            'CO₂': ['NaOH', 'KOH', 'Ca(OH)₂', 'Mg', 'C'],
            'NH₃': ['HCl', 'H₂SO₄', 'HNO₃', 'O₂'],

            // Important gases and compounds
            'H₂O₂': ['KI', 'MnO₂', 'FeCl₃'],
            'SO₂': ['O₂', 'NaOH', 'KOH', 'Ca(OH)₂'],
            'NO₂': ['H₂O', 'NaOH', 'KOH'],
            'CO': ['O₂', 'Fe₂O₃', 'CuO'],
            'H₂S': ['O₂', 'SO₂', 'FeCl₃', 'CuSO₄'],
            'N₂': ['H₂', 'O₂'],
            'P₄': ['O₂', 'Cl₂'],
            'S': ['O₂', 'Fe', 'Cu', 'Zn', 'H₂'],

            // Water
            'H₂O': ['Na', 'K', 'Ca', 'CaO', 'MgO', 'SO₃', 'P₂O₅'],
        };

        // Check each existing reactant for patterns
        existingReactants.forEach(reactant => {
            if (patterns[reactant]) {
                suggestions.push(...patterns[reactant]);
            }
        });

        // If no specific patterns found, suggest common reactants
        if (suggestions.length === 0) {
            suggestions.push('O₂', 'H₂O', 'HCl', 'NaOH', 'H₂SO₄', 'NaCl', 'AgNO₃');
        }

        // Remove duplicates and existing compounds
        return [...new Set(suggestions)].filter(s => !existingReactants.includes(s));
    }

    // Get smart product suggestions based on what's already entered
    getProductSuggestions(existingProducts) {
        const suggestions = [];

        // Common product patterns
        const patterns = {
            // If we have water, suggest salts or gases
            'H₂O': ['NaCl', 'CO₂', 'CaCl₂', 'MgCl₂'],

            // If we have CO₂, suggest water (combustion products)
            'CO₂': ['H₂O'],

            // If we have a salt, suggest water
            'NaCl': ['H₂O', 'AgCl'],
            'CaCl₂': ['H₂O', 'CO₂'],
            'MgCl₂': ['H₂O'],

            // If we have an oxide, suggest water or CO₂
            'CaO': ['H₂O', 'CO₂'],
            'MgO': ['H₂O'],
            'Fe₂O₃': ['H₂O'],
            'CuO': ['H₂O'],
            'Al₂O₃': ['H₂O'],
        };

        // Check each existing product for patterns
        existingProducts.forEach(product => {
            if (patterns[product]) {
                suggestions.push(...patterns[product]);
            }
        });

        // If no specific patterns found, suggest common products
        if (suggestions.length === 0) {
            suggestions.push('H₂O', 'CO₂', 'NaCl', 'CaCl₂', 'MgCl₂', 'H₂');
        }

        // Remove duplicates and existing compounds
        return [...new Set(suggestions)].filter(s => !existingProducts.includes(s));
    }

    // Update product suggestions based on reactants
    updateProductSuggestions() {
        const reactantsInput = document.getElementById('reactants-input');
        const productsInput = document.getElementById('products-input');

        if (!reactantsInput || !productsInput || !reactantsInput.value.trim()) {
            return;
        }

        const reactants = reactantsInput.value.split('+').map(r => r.trim()).filter(r => r);
        const suggestedProducts = this.predictProducts(reactants);

        if (suggestedProducts.length > 0 && !productsInput.value.trim()) {
            const suggestion = suggestedProducts.join(' + ');
            productsInput.value = suggestion;
            productsInput.style.color = '#9ca3af'; // Gray color
            productsInput.classList.add('has-suggestion');
        }
    }

    // Predict likely products based on reactants
    predictProducts(reactants) {
        // console.log('predictProducts called with:', reactants);
        if (reactants.length === 0) return [];

        const products = [];

        // Debug specific problematic reactions
        if (reactants.includes('H₂S') && reactants.includes('O₂')) {
            console.log('DEBUG: H₂S + O₂ detected, checking synthesis...');
        }
        if (reactants.includes('NH₃') && reactants.includes('O₂')) {
            console.log('DEBUG: NH₃ + O₂ detected, checking synthesis...');
        }
        if (reactants.includes('Mg') && reactants.includes('CO₂')) {
            console.log('DEBUG: Mg + CO₂ detected, checking reduction...');
        }

        // Handle single reactants - suggest what they commonly react with
        if (reactants.length === 1) {
            return this.getSingleReactantProducts(reactants[0]);
        }

        // Handle multiple reactants - ORDERED BY SPECIFICITY (most specific first)

        // Handle multiple reactants - ORDERED BY SPECIFICITY (most specific first)

        // 1. VERY SPECIFIC REACTIONS (exact compound matches)
        if (this.isPhotosynthesisReaction(reactants)) {
            products.push(...this.getPhotosynthesisProducts(reactants));
        } else if (this.isRespirationReaction(reactants)) {
            products.push(...this.getRespirationProducts(reactants));
        } else if (this.isFermentationReaction(reactants)) {
            products.push(...this.getFermentationProducts(reactants));
        } else if (this.isElectrolysisReaction(reactants)) {
            products.push(...this.getElectrolysisProducts(reactants));
        } else if (this.isDecompositionReaction(reactants)) {
            products.push(...this.getDecompositionProducts(reactants));
        }

        // 2. DISPLACEMENT REACTIONS (check before general compound types)
        else if (this.isDoubleDisplacement(reactants)) {
            console.log('Detected: double displacement for:', reactants);
            products.push(...this.getDoubleDisplacementProducts(reactants));
        } else if (this.isReductionReaction(reactants)) {
            console.log('Detected: reduction');
            products.push(...this.getReductionProducts(reactants));
        } else if (this.isSingleDisplacement(reactants)) {
            console.log('Single displacement detected, getting products...');
            products.push(...this.getSingleDisplacementProducts(reactants));
        } else if (this.isHalogenDisplacement(reactants)) {
            products.push(...this.getHalogenDisplacementProducts(reactants));
        } else if (this.isSynthesisReaction(reactants)) {
            console.log('Detected: synthesis reaction for:', reactants);
            products.push(...this.getSynthesisProducts(reactants));
        }

        // 2.5. SPECIFIC COMPOUND TYPE REACTIONS
        else if (this.isCarbonatAcidReaction(reactants)) {
            products.push(...this.getCarbonateAcidProducts(reactants));
        } else if (this.isAmmoniaReaction(reactants)) {
            products.push(...this.getAmmoniaProducts(reactants));
        } else if (this.isPeroxideReaction(reactants)) {
            console.log('Detected: peroxide reaction');
            products.push(...this.getPeroxideProducts(reactants));
        } else if (this.isNitrateReaction(reactants)) {
            products.push(...this.getNitrateProducts(reactants));
        } else if (this.isSulfideReaction(reactants)) {
            products.push(...this.getSulfideProducts(reactants));
        } else if (this.isPhosphateReaction(reactants)) {
            products.push(...this.getPhosphateProducts(reactants));
        }

        // 3. ORGANIC REACTIONS (specific organic patterns)
        else if (this.isHydrogenationReaction(reactants)) {
            products.push(...this.getHydrogenationProducts(reactants));
        } else if (this.isDehydrationReaction(reactants)) {
            products.push(...this.getDehydrationProducts(reactants));
            // } else if (this.isCarbideHydrolysis(reactants)) {
            //     products.push(...this.getCarbideHydrolysisProducts(reactants));
            // } else if (this.isNitrideHydrolysis(reactants)) {
            //     products.push(...this.getNitrideHydrolysisProducts(reactants));
            // } else if (this.isPolymerizationReaction(reactants)) {
            products.push(...this.getPolymerizationProducts(reactants));
        } else if (this.isOrganicReaction(reactants)) {
            products.push(...this.getOrganicProducts(reactants));
        }



        // 4. REDOX REACTIONS (specific patterns)
        else if (this.isRedoxReaction(reactants)) {
            console.log('Detected: redox reaction');
            products.push(...this.getRedoxProducts(reactants));
        }

        // 5. METAL REACTIONS (specific metal + other compound)
        else if (this.isMetalWaterReaction(reactants)) {
            products.push(...this.getMetalWaterProducts(reactants));
        } else if (this.isMetalAcidReaction(reactants)) {
            products.push(...this.getMetalAcidProducts(reactants));
        } else if (this.isMetalOxidation(reactants)) {
            products.push(...this.getMetalOxidationProducts(reactants));
        } else if (this.isComplexOxideReaction(reactants)) {
            products.push(...this.getComplexOxideProducts(reactants));
        }

        // 6. GENERAL REACTION TYPES (broader patterns)
        else if (this.isAcidBaseReaction(reactants)) {
            console.log('Detected: acid-base reaction for:', reactants);
            products.push(...this.getAcidBaseProducts(reactants));
        } else if (this.isOxideWaterReaction(reactants)) {
            products.push(...this.getOxideWaterProducts(reactants));
        } else if (this.isCombustionReaction(reactants)) {
            // console.log('Detected: combustion reaction for:', reactants);
            products.push('CO₂', 'H₂O');
        }
        // console.log('predictProducts returning:', products);
        return products;
    }

    // Get products for single reactants
    getSingleReactantProducts(reactant) {
        const singleReactantMap = {
            // Hydrocarbons
            'CH₄': ['CO₂', 'H₂O'],
            'C₂H₆': ['CO₂', 'H₂O'],
            'C₃H₈': ['CO₂', 'H₂O'],
            'C₄H₁₀': ['CO₂', 'H₂O'],
            'C₂H₄': ['CO₂', 'H₂O'],
            'C₂H₂': ['CO₂', 'H₂O'],
            'C₆H₆': ['CO₂', 'H₂O'],
            'C₂H₅OH': ['CO₂', 'H₂O'],
            'CH₃OH': ['CO₂', 'H₂O'],

            // Acids
            'HCl': ['NaCl', 'H₂O'],
            'H₂SO₄': ['Na₂SO₄', 'H₂O'],
            'HNO₃': ['NaNO₃', 'H₂O'],
            'CH₃COOH': ['CH₃COONa', 'H₂O'],
            'H₃PO₄': ['Na₃PO₄', 'H₂O'],

            // Bases
            'NaOH': ['NaCl', 'H₂O'],
            'KOH': ['KCl', 'H₂O'],
            'Ca(OH)₂': ['CaCl₂', 'H₂O'],
            'Mg(OH)₂': ['MgCl₂', 'H₂O'],
            'Ba(OH)₂': ['BaCl₂', 'H₂O'],

            // Metals
            'Fe': ['Fe₂O₃'],
            'Cu': ['CuO'],
            'Zn': ['ZnO'],
            'Al': ['Al₂O₃'],
            'Mg': ['MgO'],
            'Ca': ['CaO'],
            'Na': ['Na₂O'],
            'K': ['K₂O'],

            // Elements
            'H₂': ['H₂O'],
            'O₂': ['H₂O'],
            'N₂': ['NH₃'],
            'Cl₂': ['NaCl'],

            // Compounds that decompose
            'H₂O₂': ['H₂O', 'O₂'],
            'KClO₃': ['KCl', 'O₂'],
            'CaCO₃': ['CaO', 'CO₂'],
            'H₂CO₃': ['H₂O', 'CO₂'],
            'NH₄Cl': ['NH₃', 'HCl'],

            // Metal oxides with reducing agents
            'Fe₂O₃': ['Fe', 'CO₂'], // with CO
            'CuO': ['Cu', 'H₂O'], // with H₂
            'ZnO': ['Zn', 'CO'], // with C
            'PbO': ['Pb', 'CO₂'], // with CO
            'MnO₂': ['Mn', 'CO₂'], // with C

            // Steam reactions (high temperature metal + water)
            'Fe': ['Fe₃O₄', 'H₂'], // with H₂O at high temp
            'Zn': ['ZnO', 'H₂'], // with H₂O at high temp
            'Al': ['Al₂O₃', 'H₂'], // with H₂O at high temp
            'Cr': ['Cr₂O₃', 'H₂'], // with H₂O at high temp
            'Ni': ['NiO', 'H₂'], // with H₂O at high temp

            // Complex compounds
            'KMnO₄': ['KCl', 'MnCl₂', 'Cl₂', 'H₂O'], // with HCl
            'AgNO₃': ['AgCl', 'NaNO₃'], // with NaCl
            'Pb(NO₃)₂': ['PbI₂', 'KNO₃'], // with KI
            'FeCl₃': ['FeCl₂', 'KCl', 'I₂'], // with KI

            // Biological reactions
            'C₆H₁₂O₆': ['CO₂', 'H₂O'], // respiration with O₂
            'CO₂': ['C₆H₁₂O₆', 'O₂'], // photosynthesis with H₂O

            // Organic reactions
            'C₂H₄': ['C₂H₆'], // hydrogenation with H₂
            'C₂H₂': ['C₂H₄'], // hydrogenation with H₂
            'C₂H₅OH': ['C₂H₄', 'H₂O'], // dehydration
            'CH₃OH': ['HCHO', 'H₂O'], // oxidation

            // Electrolysis
            'H₂O': ['H₂', 'O₂'], // electrolysis
            'NaCl': ['Na', 'Cl₂'], // electrolysis
            'CuSO₄': ['Cu', 'O₂', 'H₂SO₄'], // electrolysis
            'Al₂O₃': ['Al', 'O₂'], // electrolysis
        };

        return singleReactantMap[reactant] || [];
    }

    // Check for synthesis reactions (A + B → AB)
    isSynthesisReaction(reactants) {
        if (reactants.length !== 2) return false;

        const synthesisPairs = [
            ['H₂', 'O₂'], ['H₂', 'Cl₂'], ['H₂', 'Br₂'], ['H₂', 'I₂'],
            ['Na', 'Cl₂'], ['K', 'Cl₂'], ['Ca', 'Cl₂'], ['Mg', 'Cl₂'],
            ['N₂', 'H₂'], ['N₂', 'O₂'], ['S', 'O₂'], ['P', 'O₂'], ['P₄', 'O₂'],
            ['C', 'O₂'], ['CO', 'O₂'], ['H₂S', 'O₂'], ['NH₃', 'O₂'], ['SO₂', 'O₂'],
            ['Fe', 'S'], ['Cu', 'S'], ['Zn', 'S'], ['Pb', 'S']
        ];

        return synthesisPairs.some(pair =>
            (reactants.includes(pair[0]) && reactants.includes(pair[1])) ||
            (reactants.includes(pair[1]) && reactants.includes(pair[0]))
        );
    }

    getSynthesisProducts(reactants) {
        const synthesisMap = {
            'H₂+O₂': ['H₂O'],
            'Cl₂+H₂': ['HCl'],
            'Br₂+H₂': ['HBr'],
            'H₂+I₂': ['HI'],
            'Cl₂+Na': ['NaCl'],
            'Cl₂+K': ['KCl'],
            'Ca+Cl₂': ['CaCl₂'],
            'Cl₂+Mg': ['MgCl₂'],
            'H₂+N₂': ['NH₃'],
            'N₂+O₂': ['NO'],
            'O₂+S': ['SO₂'],
            'O₂+P': ['P₂O₅'],
            'O₂+P₄': ['P₂O₅'],
            'C+O₂': ['CO₂'],
            'CO+O₂': ['CO₂'],
            'H₂S+O₂': ['SO₂', 'H₂O'],
            'NH₃+O₂': ['NO', 'H₂O'],
            'O₂+SO₂': ['SO₃'],
            'Fe+S': ['FeS'],
            'Cu+S': ['CuS'],
            'Pb+S': ['PbS'],
            'S+Zn': ['ZnS']
        };

        const key = reactants.sort().join('+');
        console.log('Synthesis key generated:', key, 'Found in map:', synthesisMap[key]);
        return synthesisMap[key] || [];
    }

    // Check for decomposition reactions
    isDecompositionReaction(reactants) {
        if (reactants.length !== 1) return false;

        const decomposableCompounds = [
            'H₂O₂', 'KClO₃', 'NaClO₃', 'CaCO₃', 'MgCO₃', 'Na₂CO₃',
            'H₂CO₃', 'NH₄Cl', 'NH₄NO₃', 'CuCO₃', 'ZnCO₃', 'PbCO₃',
            'Ag₂O', 'HgO', 'Cu(OH)₂', 'Fe(OH)₃', 'Al(OH)₃'
        ];

        return decomposableCompounds.includes(reactants[0]);
    }

    getDecompositionProducts(reactants) {
        const decompositionMap = {
            'H₂O₂': ['H₂O', 'O₂'],
            'KClO₃': ['KCl', 'O₂'],
            'NaClO₃': ['NaCl', 'O₂'],
            'CaCO₃': ['CaO', 'CO₂'],
            'MgCO₃': ['MgO', 'CO₂'],
            'Na₂CO₃': ['Na₂O', 'CO₂'],
            'H₂CO₃': ['H₂O', 'CO₂'],
            'NH₄Cl': ['NH₃', 'HCl'],
            'NH₄NO₃': ['N₂O', 'H₂O'],
            'CuCO₃': ['CuO', 'CO₂'],
            'ZnCO₃': ['ZnO', 'CO₂'],
            'PbCO₃': ['PbO', 'CO₂'],
            'Ag₂O': ['Ag', 'O₂'],
            'HgO': ['Hg', 'O₂'],
            'Cu(OH)₂': ['CuO', 'H₂O'],
            'Fe(OH)₃': ['Fe₂O₃', 'H₂O'],
            'Al(OH)₃': ['Al₂O₃', 'H₂O']
        };

        return decompositionMap[reactants[0]] || [];
    }

    // Check for single displacement reactions
    isSingleDisplacement(reactants) {
        if (reactants.length !== 2) return false;

        const metals = ['Fe', 'Cu', 'Zn', 'Al', 'Mg', 'Ca', 'Na', 'K', 'Li', 'Ag', 'Pb'];
        const metalSalts = ['CuSO₄', 'ZnSO₄', 'FeSO₄', 'AgNO₃', 'Pb(NO₃)₂', 'CuCl₂', 'FeCl₃', 'AlCl₃', 'Fe₂O₃'];

        const hasMetal = reactants.some(r => metals.includes(r));
        const hasMetalSalt = reactants.some(r => metalSalts.includes(r));

        console.log('isSingleDisplacement check:', reactants, 'hasMetal:', hasMetal, 'hasMetalSalt:', hasMetalSalt);

        // Single displacement requires a free metal + a metal salt/oxide (not water, acids, etc.)
        return hasMetal && hasMetalSalt;
    }

    getSingleDisplacementProducts(reactants) {
        console.log('getSingleDisplacementProducts called with:', reactants);
        // Common displacement reactions based on activity series
        if (reactants.includes('Zn') && reactants.includes('CuSO₄')) {
            console.log('Found Zn + CuSO₄ reaction');
            return ['ZnSO₄', 'Cu'];
        }
        if (reactants.includes('Fe') && reactants.includes('CuSO₄')) {
            console.log('Found Fe + CuSO₄ reaction');
            return ['FeSO₄', 'Cu'];
        }
        if (reactants.includes('Al') && reactants.includes('Fe₂O₃')) {
            console.log('Found Al + Fe₂O₃ reaction');
            return ['Al₂O₃', 'Fe'];
        }
        if (reactants.includes('Cu') && reactants.includes('AgNO₃')) {
            console.log('Found Cu + AgNO₃ reaction');
            return ['Cu(NO₃)₂', 'Ag'];
        }
        if (reactants.includes('Zn') && reactants.includes('AgNO₃')) {
            console.log('Found Zn + AgNO₃ reaction');
            return ['Zn(NO₃)₂', 'Ag'];
        }
        if (reactants.includes('Fe') && reactants.includes('AgNO₃')) {
            console.log('Found Fe + AgNO₃ reaction');
            return ['Fe(NO₃)₃', 'Ag'];
        }
        console.log('No single displacement reaction found for:', reactants);
        return [];
    }

    // Metal + Acid reactions
    isMetalAcidReaction(reactants) {
        const metals = ['Fe', 'Cu', 'Zn', 'Al', 'Mg', 'Ca', 'Na', 'K'];
        const acids = ['HCl', 'H₂SO₄', 'HNO₃'];

        return reactants.some(r => metals.includes(r)) &&
            reactants.some(r => acids.includes(r));
    }

    getMetalAcidProducts(reactants) {
        // Metal + Acid → Salt + H₂
        if (reactants.includes('Zn') && reactants.includes('HCl')) {
            return ['ZnCl₂', 'H₂'];
        }
        if (reactants.includes('Mg') && reactants.includes('HCl')) {
            return ['MgCl₂', 'H₂'];
        }
        if (reactants.includes('Fe') && reactants.includes('HCl')) {
            return ['FeCl₂', 'H₂'];
        }
        if (reactants.includes('Al') && reactants.includes('HCl')) {
            return ['AlCl₃', 'H₂'];
        }
        return ['H₂']; // Generic hydrogen gas production
    }

    // Metal + Water reactions (expanded to include steam reactions)
    isMetalWaterReaction(reactants) {
        if (reactants.length !== 2) return false;

        const metals = ['Na', 'K', 'Li', 'Ca', 'Mg', 'Fe', 'Zn', 'Al', 'Cr', 'Ni', 'Co', 'Mn', 'Sn', 'Pb'];
        const hasWater = reactants.includes('H₂O');
        const hasMetal = reactants.some(r => metals.includes(r));

        return hasWater && hasMetal;
    }

    getMetalWaterProducts(reactants) {
        // Active metals at room temperature
        if (reactants.includes('Na')) return ['NaOH', 'H₂'];
        if (reactants.includes('K')) return ['KOH', 'H₂'];
        if (reactants.includes('Li')) return ['LiOH', 'H₂'];
        if (reactants.includes('Ca')) return ['Ca(OH)₂', 'H₂'];
        if (reactants.includes('Mg')) return ['Mg(OH)₂', 'H₂'];

        // Steam reactions (high temperature)
        if (reactants.includes('Fe')) return ['Fe₃O₄', 'H₂'];
        if (reactants.includes('Zn')) return ['ZnO', 'H₂'];
        if (reactants.includes('Al')) return ['Al₂O₃', 'H₂'];
        if (reactants.includes('Cr')) return ['Cr₂O₃', 'H₂'];
        if (reactants.includes('Ni')) return ['NiO', 'H₂'];
        if (reactants.includes('Co')) return ['CoO', 'H₂'];
        if (reactants.includes('Mn')) return ['MnO', 'H₂'];
        if (reactants.includes('Sn')) return ['SnO₂', 'H₂'];
        if (reactants.includes('Pb')) return ['PbO', 'H₂'];

        return ['H₂'];
    }

    // Oxide + Water reactions
    isOxideWaterReaction(reactants) {
        const oxides = ['CaO', 'MgO', 'Na₂O', 'K₂O', 'SO₃', 'CO₂', 'P₂O₅'];
        return reactants.includes('H₂O') &&
            reactants.some(r => oxides.includes(r));
    }

    getOxideWaterProducts(reactants) {
        if (reactants.includes('CaO')) return ['Ca(OH)₂'];
        if (reactants.includes('MgO')) return ['Mg(OH)₂'];
        if (reactants.includes('Na₂O')) return ['NaOH'];
        if (reactants.includes('K₂O')) return ['KOH'];
        if (reactants.includes('SO₃')) return ['H₂SO₄'];
        if (reactants.includes('CO₂')) return ['H₂CO₃'];
        if (reactants.includes('P₂O₅')) return ['H₃PO₄'];
        return [];
    }

    // Carbonate + Acid reactions
    isCarbonatAcidReaction(reactants) {
        const carbonates = ['CaCO₃', 'MgCO₃', 'Na₂CO₃', 'K₂CO₃', 'ZnCO₃', 'CuCO₃'];
        const acids = ['HCl', 'H₂SO₄', 'HNO₃', 'CH₃COOH'];

        return reactants.some(r => carbonates.includes(r)) &&
            reactants.some(r => acids.includes(r));
    }

    getCarbonateAcidProducts(reactants) {
        // Carbonate + Acid → Salt + CO₂ + H₂O
        if (reactants.includes('CaCO₃') && reactants.includes('HCl')) {
            return ['CaCl₂', 'CO₂', 'H₂O'];
        }
        if (reactants.includes('Na₂CO₃') && reactants.includes('HCl')) {
            return ['NaCl', 'CO₂', 'H₂O'];
        }
        return ['CO₂', 'H₂O']; // Generic carbonate reaction
    }

    // Reduction reactions (Metal oxide + reducing agent)
    isReductionReaction(reactants) {
        const metalOxides = ['Fe₂O₃', 'CuO', 'ZnO', 'Al₂O₃', 'MgO', 'CaO', 'PbO', 'SnO₂', 'Cr₂O₃', 'MnO₂', 'CO₂'];
        const reducingAgents = ['CO', 'H₂', 'C', 'Al', 'Mg', 'Zn', 'Fe'];

        return reactants.some(r => metalOxides.includes(r)) &&
            reactants.some(r => reducingAgents.includes(r));
    }

    getReductionProducts(reactants) {
        // Fe₂O₃ + CO → Fe + CO₂
        if (reactants.includes('Fe₂O₃') && reactants.includes('CO')) {
            return ['Fe', 'CO₂'];
        }
        if (reactants.includes('CuO') && reactants.includes('H₂')) {
            return ['Cu', 'H₂O'];
        }
        if (reactants.includes('ZnO') && reactants.includes('C')) {
            return ['Zn', 'CO'];
        }
        if (reactants.includes('Al₂O₃') && reactants.includes('C')) {
            return ['Al', 'CO'];
        }
        if (reactants.includes('PbO') && reactants.includes('CO')) {
            return ['Pb', 'CO₂'];
        }
        if (reactants.includes('CuO') && reactants.includes('Al')) {
            return ['Cu', 'Al₂O₃'];
        }
        if (reactants.includes('CO₂') && reactants.includes('Mg')) {
            return ['MgO', 'CO'];
        }
        if (reactants.includes('CuO') && reactants.includes('CO')) {
            return ['Cu', 'CO₂'];
        }
        if (reactants.includes('Fe₂O₃') && reactants.includes('H₂')) {
            return ['Fe', 'H₂O'];
        }
        if (reactants.includes('Fe₂O₃') && reactants.includes('C')) {
            return ['Fe', 'CO'];
        }
        if (reactants.includes('MnO₂') && reactants.includes('C')) {
            return ['Mn', 'CO₂'];
        }
        if (reactants.includes('Cr₂O₃') && reactants.includes('Al')) {
            return ['Cr', 'Al₂O₃'];
        }
        return [];
    }

    // Redox reactions (more complex electron transfer)
    isRedoxReaction(reactants) {
        const redoxPairs = [
            ['KMnO₄', 'HCl'], ['K₂Cr₂O₇', 'HCl'], ['H₂O₂', 'KI'],
            ['Br₂', 'KI'], ['Cl₂', 'KBr'], ['I₂', 'Na₂S₂O₃'],
            ['FeCl₃', 'KI'], ['CuSO₄', 'KI'], ['AgNO₃', 'Cu']
        ];

        return redoxPairs.some(pair =>
            (reactants.includes(pair[0]) && reactants.includes(pair[1])) ||
            (reactants.includes(pair[1]) && reactants.includes(pair[0]))
        );
    }

    getRedoxProducts(reactants) {
        if (reactants.includes('KMnO₄') && reactants.includes('HCl')) {
            return ['KCl', 'MnCl₂', 'Cl₂', 'H₂O'];
        }
        if (reactants.includes('K₂Cr₂O₇') && reactants.includes('HCl')) {
            return ['KCl', 'CrCl₃', 'Cl₂', 'H₂O'];
        }
        if (reactants.includes('H₂O₂') && reactants.includes('KI')) {
            return ['KOH', 'I₂', 'H₂O'];
        }
        if (reactants.includes('Br₂') && reactants.includes('KI')) {
            return ['KBr', 'I₂'];
        }
        if (reactants.includes('Cl₂') && reactants.includes('KBr')) {
            return ['KCl', 'Br₂'];
        }
        if (reactants.includes('I₂') && reactants.includes('Na₂S₂O₃')) {
            return ['NaI', 'Na₂S₄O₆'];
        }
        if (reactants.includes('FeCl₃') && reactants.includes('KI')) {
            return ['FeCl₂', 'KCl', 'I₂'];
        }
        if (reactants.includes('CuSO₄') && reactants.includes('KI')) {
            return ['CuI', 'K₂SO₄', 'I₂'];
        }
        if (reactants.includes('AgNO₃') && reactants.includes('Cu')) {
            return ['Ag', 'Cu(NO₃)₂'];
        }
        return [];
    }

    // Halogen displacement reactions
    isHalogenDisplacement(reactants) {
        const halogens = ['F₂', 'Cl₂', 'Br₂', 'I₂'];
        const halides = ['NaF', 'NaCl', 'NaBr', 'NaI', 'KF', 'KCl', 'KBr', 'KI'];

        return reactants.some(r => halogens.includes(r)) &&
            reactants.some(r => halides.includes(r));
    }

    getHalogenDisplacementProducts(reactants) {
        // Sodium halide reactions
        if (reactants.includes('Cl₂') && reactants.includes('NaBr')) {
            return ['NaCl', 'Br₂'];
        }
        if (reactants.includes('Cl₂') && reactants.includes('NaI')) {
            return ['NaCl', 'I₂'];
        }
        if (reactants.includes('Br₂') && reactants.includes('NaI')) {
            return ['NaBr', 'I₂'];
        }
        if (reactants.includes('F₂') && reactants.includes('NaCl')) {
            return ['NaF', 'Cl₂'];
        }

        // Potassium halide reactions
        if (reactants.includes('Cl₂') && reactants.includes('KBr')) {
            return ['KCl', 'Br₂'];
        }
        if (reactants.includes('Cl₂') && reactants.includes('KI')) {
            return ['KCl', 'I₂'];
        }
        if (reactants.includes('Br₂') && reactants.includes('KI')) {
            return ['KBr', 'I₂'];
        }
        if (reactants.includes('F₂') && reactants.includes('KCl')) {
            return ['KF', 'Cl₂'];
        }
        if (reactants.includes('F₂') && reactants.includes('KBr')) {
            return ['KF', 'Br₂'];
        }
        if (reactants.includes('F₂') && reactants.includes('KI')) {
            return ['KF', 'I₂'];
        }

        return [];
    }

    // Ammonia-related reactions
    isAmmoniaReaction(reactants) {
        return reactants.includes('NH₃') || reactants.includes('NH₄Cl') ||
            reactants.includes('(NH₄)₂SO₄') || reactants.includes('NH₄NO₃');
    }

    getAmmoniaProducts(reactants) {
        if (reactants.includes('NH₃') && reactants.includes('HCl')) {
            return ['NH₄Cl'];
        }
        if (reactants.includes('NH₃') && reactants.includes('H₂SO₄')) {
            return ['(NH₄)₂SO₄'];
        }
        if (reactants.includes('NH₃') && reactants.includes('HNO₃')) {
            return ['NH₄NO₃'];
        }
        if (reactants.includes('NH₄Cl') && reactants.includes('NaOH')) {
            return ['NaCl', 'NH₃', 'H₂O'];
        }
        if (reactants.includes('(NH₄)₂SO₄') && reactants.includes('Ca(OH)₂')) {
            return ['CaSO₄', 'NH₃', 'H₂O'];
        }
        return [];
    }

    // Peroxide reactions
    isPeroxideReaction(reactants) {
        const peroxides = ['H₂O₂', 'Na₂O₂', 'BaO₂'];
        return reactants.some(r => peroxides.includes(r));
    }

    getPeroxideProducts(reactants) {
        if (reactants.includes('H₂O₂') && reactants.includes('KI')) {
            return ['KOH', 'I₂'];
        }
        if (reactants.includes('H₂O₂') && reactants.includes('MnO₂')) {
            return ['H₂O', 'O₂', 'MnO₂']; // MnO₂ is catalyst
        }
        if (reactants.includes('H₂O₂') && reactants.includes('FeCl₃')) {
            return ['FeCl₂', 'HCl', 'H₂O', 'O₂'];
        }
        if (reactants.includes('Na₂O₂') && reactants.includes('H₂O')) {
            return ['NaOH', 'O₂'];
        }
        if (reactants.includes('BaO₂') && reactants.includes('H₂SO₄')) {
            return ['BaSO₄', 'H₂O₂'];
        }
        return [];
    }

    // Nitrate reactions
    isNitrateReaction(reactants) {
        const nitrates = ['AgNO₃', 'Pb(NO₃)₂', 'Cu(NO₃)₂', 'Fe(NO₃)₃', 'KNO₃', 'NaNO₃'];
        return reactants.some(r => nitrates.includes(r));
    }

    getNitrateProducts(reactants) {
        if (reactants.includes('AgNO₃') && reactants.includes('NaCl')) {
            return ['AgCl', 'NaNO₃'];
        }
        if (reactants.includes('Pb(NO₃)₂') && reactants.includes('KI')) {
            return ['PbI₂', 'KNO₃'];
        }
        if (reactants.includes('AgNO₃') && reactants.includes('KBr')) {
            return ['AgBr', 'KNO₃'];
        }
        if (reactants.includes('Cu(NO₃)₂') && reactants.includes('NaOH')) {
            return ['Cu(OH)₂', 'NaNO₃'];
        }
        if (reactants.includes('Fe(NO₃)₃') && reactants.includes('KOH')) {
            return ['Fe(OH)₃', 'KNO₃'];
        }
        return [];
    }

    // Sulfide reactions
    isSulfideReaction(reactants) {
        const sulfides = ['H₂S', 'Na₂S', 'FeS', 'CuS', 'ZnS', 'PbS'];
        return reactants.some(r => sulfides.includes(r));
    }

    getSulfideProducts(reactants) {
        if (reactants.includes('H₂S') && reactants.includes('SO₂')) {
            return ['S', 'H₂O'];
        }
        if (reactants.includes('Na₂S') && reactants.includes('HCl')) {
            return ['NaCl', 'H₂S'];
        }
        if (reactants.includes('FeS') && reactants.includes('HCl')) {
            return ['FeCl₂', 'H₂S'];
        }
        if (reactants.includes('CuS') && reactants.includes('HNO₃')) {
            return ['Cu(NO₃)₂', 'S', 'NO', 'H₂O'];
        }
        return [];
    }

    // Phosphate reactions
    isPhosphateReaction(reactants) {
        const phosphates = ['H₃PO₄', 'Na₃PO₄', 'Ca₃(PO₄)₂', 'K₃PO₄'];
        return reactants.some(r => phosphates.includes(r));
    }

    getPhosphateProducts(reactants) {
        if (reactants.includes('H₃PO₄') && reactants.includes('NaOH')) {
            return ['Na₃PO₄', 'H₂O'];
        }
        if (reactants.includes('Ca₃(PO₄)₂') && reactants.includes('H₂SO₄')) {
            return ['CaSO₄', 'H₃PO₄'];
        }
        if (reactants.includes('Na₃PO₄') && reactants.includes('AgNO₃')) {
            return ['Ag₃PO₄', 'NaNO₃'];
        }
        return [];
    }

    // Complex oxide formation reactions
    isComplexOxideReaction(reactants) {
        const metals = ['Fe', 'Cr', 'Mn', 'Co', 'Ni', 'Cu', 'Pb', 'Sn'];
        const oxidizers = ['O₂', 'H₂O', 'CO₂', 'NO₂', 'SO₂'];

        return reactants.some(r => metals.includes(r)) &&
            reactants.some(r => oxidizers.includes(r));
    }

    getComplexOxideProducts(reactants) {
        // Mixed oxidation states
        if (reactants.includes('Fe') && reactants.includes('O₂')) {
            return ['Fe₃O₄']; // Magnetite (mixed Fe²⁺/Fe³⁺)
        }
        if (reactants.includes('Pb') && reactants.includes('O₂')) {
            return ['Pb₃O₄']; // Red lead
        }
        if (reactants.includes('Mn') && reactants.includes('O₂')) {
            return ['Mn₃O₄']; // Hausmannite
        }
        if (reactants.includes('Co') && reactants.includes('O₂')) {
            return ['Co₃O₄']; // Cobalt oxide
        }
        return [];
    }

    // Organic reactions
    isOrganicReaction(reactants) {
        const organics = ['C₂H₄', 'C₂H₂', 'CH₃OH', 'C₂H₅OH', 'CH₃COOH', 'C₆H₆', 'C₆H₁₂O₆'];
        const reagents = ['H₂', 'Br₂', 'Cl₂', 'H₂O', 'O₂', 'KMnO₄', 'K₂Cr₂O₇'];

        return reactants.some(r => organics.includes(r)) &&
            reactants.some(r => reagents.includes(r));
    }

    getOrganicProducts(reactants) {
        // Addition reactions
        if (reactants.includes('C₂H₄') && reactants.includes('H₂')) {
            return ['C₂H₆'];
        }
        if (reactants.includes('C₂H₄') && reactants.includes('Br₂')) {
            return ['C₂H₄Br₂'];
        }
        if (reactants.includes('C₂H₄') && reactants.includes('H₂O')) {
            return ['C₂H₅OH'];
        }
        if (reactants.includes('C₂H₂') && reactants.includes('H₂')) {
            return ['C₂H₄'];
        }
        if (reactants.includes('C₂H₂') && reactants.includes('H₂O')) {
            return ['CH₃CHO'];
        }

        // Substitution reactions
        if (reactants.includes('C₆H₆') && reactants.includes('Br₂')) {
            return ['C₆H₅Br', 'HBr'];
        }
        if (reactants.includes('C₆H₆') && reactants.includes('Cl₂')) {
            return ['C₆H₅Cl', 'HCl'];
        }

        // Oxidation reactions
        if (reactants.includes('C₂H₅OH') && reactants.includes('O₂')) {
            return ['CH₃COOH', 'H₂O'];
        }
        if (reactants.includes('C₂H₅OH') && reactants.includes('K₂Cr₂O₇')) {
            return ['CH₃COOH', 'Cr₂O₃', 'K₂O', 'H₂O'];
        }
        if (reactants.includes('CH₃OH') && reactants.includes('KMnO₄')) {
            return ['HCHO', 'MnO₂', 'KOH', 'H₂O'];
        }
        if (reactants.includes('CH₃OH') && reactants.includes('O₂')) {
            return ['HCHO', 'H₂O'];
        }

        return [];
    }

    // Electrolysis reactions
    isElectrolysisReaction(reactants) {
        const electrolytes = ['H₂O', 'NaCl', 'CuSO₄', 'Al₂O₃', 'KBr', 'AgNO₃'];
        return reactants.length === 1 && electrolytes.includes(reactants[0]);
    }

    getElectrolysisProducts(reactants) {
        if (reactants.includes('H₂O')) {
            return ['H₂', 'O₂'];
        }
        if (reactants.includes('NaCl')) {
            return ['Na', 'Cl₂'];
        }
        if (reactants.includes('CuSO₄')) {
            return ['Cu', 'O₂', 'H₂SO₄'];
        }
        if (reactants.includes('Al₂O₃')) {
            return ['Al', 'O₂'];
        }
        if (reactants.includes('KBr')) {
            return ['K', 'Br₂'];
        }
        return [];
    }

    // Photosynthesis reaction
    isPhotosynthesisReaction(reactants) {
        return reactants.includes('CO₂') && reactants.includes('H₂O');
    }

    getPhotosynthesisProducts(reactants) {
        if (reactants.includes('CO₂') && reactants.includes('H₂O')) {
            return ['C₆H₁₂O₆', 'O₂'];
        }
        return [];
    }

    // Cellular respiration
    isRespirationReaction(reactants) {
        return reactants.includes('C₆H₁₂O₆') && reactants.includes('O₂');
    }

    getRespirationProducts(reactants) {
        if (reactants.includes('C₆H₁₂O₆') && reactants.includes('O₂')) {
            return ['CO₂', 'H₂O'];
        }
        return [];
    }

    // Fermentation reactions
    isFermentationReaction(reactants) {
        return reactants.includes('C₆H₁₂O₆') && !reactants.includes('O₂');
    }

    getFermentationProducts(reactants) {
        if (reactants.includes('C₆H₁₂O₆')) {
            return ['C₂H₅OH', 'CO₂']; // Alcoholic fermentation
        }
        return [];
    }

    // Polymerization reactions
    isPolymerizationReaction(reactants) {
        const monomers = ['C₂H₄', 'C₃H₆', 'C₄H₆', 'C₈H₈'];
        return reactants.some(r => monomers.includes(r));
    }

    getPolymerizationProducts(reactants) {
        if (reactants.includes('C₂H₄')) {
            return ['(C₂H₄)ₙ']; // Polyethylene
        }
        if (reactants.includes('C₃H₆')) {
            return ['(C₃H₆)ₙ']; // Polypropylene
        }
        if (reactants.includes('C₈H₈')) {
            return ['(C₈H₈)ₙ']; // Polystyrene
        }
        return [];
    }

    // Hydrogenation reactions
    isHydrogenationReaction(reactants) {
        const unsaturated = ['C₂H₄', 'C₂H₂', 'C₃H₆', 'C₄H₆'];
        return reactants.includes('H₂') &&
            reactants.some(r => unsaturated.includes(r));
    }

    getHydrogenationProducts(reactants) {
        if (reactants.includes('C₂H₄') && reactants.includes('H₂')) {
            return ['C₂H₆'];
        }
        if (reactants.includes('C₂H₂') && reactants.includes('H₂')) {
            return ['C₂H₄'];
        }
        if (reactants.includes('C₃H₆') && reactants.includes('H₂')) {
            return ['C₃H₈'];
        }
        return [];
    }

    // Dehydration reactions
    isDehydrationReaction(reactants) {
        const alcohols = ['C₂H₅OH', 'CH₃OH', 'C₃H₇OH'];
        const dehydrators = ['H₂SO₄', 'Al₂O₃'];

        return reactants.some(r => alcohols.includes(r)) &&
            reactants.some(r => dehydrators.includes(r));
    }

    getDehydrationProducts(reactants) {
        if (reactants.includes('C₂H₅OH') && reactants.includes('H₂SO₄')) {
            return ['C₂H₄', 'H₂O'];
        }
        if (reactants.includes('CH₃OH') && reactants.includes('Al₂O₃')) {
            return ['CH₂O', 'H₂O'];
        }
        return [];
    }

    // Helper functions for single reactant predictions
    isHydrocarbon(compound) {
        const hydrocarbons = ['CH₄', 'C₂H₆', 'C₃H₈', 'C₄H₁₀', 'C₂H₄', 'C₂H₂', 'C₆H₆', 'C₂H₅OH', 'CH₃OH'];
        return hydrocarbons.includes(compound);
    }

    isAcid(compound) {
        const acids = ['HCl', 'H₂SO₄', 'HNO₃', 'CH₃COOH', 'H₃PO₄', 'HF', 'HBr', 'HI'];
        return acids.includes(compound);
    }

    isBase(compound) {
        const bases = ['NaOH', 'KOH', 'Ca(OH)₂', 'Mg(OH)₂', 'Ba(OH)₂', 'Al(OH)₃', 'Fe(OH)₂', 'Fe(OH)₃'];
        return bases.includes(compound);
    }

    isMetal(compound) {
        const metals = ['Fe', 'Cu', 'Zn', 'Al', 'Mg', 'Ca', 'Na', 'K', 'Li', 'Ag', 'Pb', 'Sn', 'Ni', 'Co', 'Mn', 'Cr'];
        return metals.includes(compound);
    }

    getMetalOxide(metal) {
        const oxides = {
            'Fe': 'Fe₂O₃',
            'Cu': 'CuO',
            'Zn': 'ZnO',
            'Al': 'Al₂O₃',
            'Mg': 'MgO',
            'Ca': 'CaO',
            'Na': 'Na₂O',
            'K': 'K₂O',
            'Li': 'Li₂O',
            'Ag': 'Ag₂O',
            'Pb': 'PbO',
            'Sn': 'SnO₂',
            'Ni': 'NiO',
            'Co': 'CoO',
            'Mn': 'MnO',
            'Cr': 'Cr₂O₃'
        };
        return oxides[metal] || metal + 'O';
    }

    // Show smart suggestions based on content in the other field
    showSmartSuggestionsBasedOnOtherField(targetType, otherFieldValue) {
        console.log('showSmartSuggestionsBasedOnOtherField called:', targetType, otherFieldValue);

        const compounds = otherFieldValue.split('+').map(c => c.trim()).filter(c => c);
        console.log('Parsed compounds:', compounds);

        let suggestions = [];

        if (targetType === 'products') {
            // Show products based on reactants
            console.log('Getting products for reactants:', compounds);
            suggestions = this.predictProducts(compounds);
            console.log('Product suggestions:', suggestions);
        } else {
            // Show reactants based on products
            console.log('Getting reactants for products:', compounds);
            suggestions = this.predictReactants(compounds);
            console.log('Reactant suggestions:', suggestions);
        }

        if (suggestions.length === 0) {
            console.log('No suggestions found, showing regular dropdown');
            // Fallback to regular dropdown
            this.showDropdown(targetType);
            this.filterDropdown(targetType, '');
            return;
        }

        console.log('Found suggestions, creating dropdown with:', suggestions);

        const listElement = document.getElementById(`${targetType}-list`);
        if (!listElement) return;

        console.log('Showing smart suggestions based on other field:', suggestions);
        listElement.innerHTML = '';

        // Add header for smart suggestions
        const header = document.createElement('div');
        header.className = 'px-3 py-2 text-xs font-semibold text-gray-500 uppercase tracking-wide suggestion-header';
        header.textContent = '💡 Recommended for this equation';
        listElement.appendChild(header);

        // Create complete equation option instead of individual compounds
        if (suggestions.length > 0) {
            const completeEquation = suggestions.join(' + ');
            const option = this.createCompleteEquationOption(completeEquation, targetType, suggestions);
            option.classList.add('smart-suggestion');
            listElement.appendChild(option);
        }

        // Add separator
        const separator = document.createElement('div');
        separator.className = 'border-t border-gray-200 my-2';
        listElement.appendChild(separator);

        // Add "All Compounds" header
        const allHeader = document.createElement('div');
        allHeader.className = 'px-3 py-2 text-xs font-semibold text-gray-500 uppercase tracking-wide bg-gray-50';
        allHeader.textContent = 'All Compounds';
        listElement.appendChild(allHeader);

        // Group remaining compounds by category
        const categories = {};
        this.compounds.forEach(compound => {
            if (!categories[compound.category]) {
                categories[compound.category] = [];
            }
            categories[compound.category].push(compound);
        });

        // Add all compounds organized by category (excluding already suggested ones)
        Object.keys(categories).sort().forEach(category => {
            // Add category header
            const categoryHeader = document.createElement('div');
            categoryHeader.className = 'px-3 py-1 text-xs font-medium text-gray-600 bg-gray-50';
            categoryHeader.textContent = category;
            listElement.appendChild(categoryHeader);

            // Add compounds in this category (skip if already in suggestions)
            categories[category].forEach(compound => {
                if (suggestions.includes(compound.formula)) return;

                const option = this.createCompoundOption(compound, targetType);
                listElement.appendChild(option);
            });
        });

        this.showDropdown(targetType);
    }

    // Check if it's an acid-base reaction
    isAcidBaseReaction(reactants) {
        if (reactants.length !== 2) return false;

        const acids = ['HCl', 'H₂SO₄', 'HNO₃', 'CH₃COOH', 'H₃PO₄', 'HF', 'HBr', 'HI', 'SO₂', 'NO₂'];
        const bases = ['NaOH', 'KOH', 'Ca(OH)₂', 'Mg(OH)₂', 'Ba(OH)₂', 'Al(OH)₃', 'LiOH', 'H₂O'];
        const salts = ['MgSO₄', 'ZnSO₄', 'FeSO₄', 'AlCl₃', 'MgCl₂', 'BaCl₂'];

        const hasAcid = reactants.some(r => acids.includes(r));
        const hasBase = reactants.some(r => bases.includes(r));
        const hasSalt = reactants.some(r => salts.includes(r));

        // Can be acid + base, salt + base, or acid + salt
        return (hasAcid && hasBase) || (hasSalt && hasBase) || (hasAcid && hasSalt);
    }

    // Get acid-base reaction products
    getAcidBaseProducts(reactants) {
        // Common acid-base reactions
        if (reactants.includes('HCl') && reactants.includes('NaOH')) {
            return ['NaCl', 'H₂O'];
        }
        if (reactants.includes('HCl') && reactants.includes('KOH')) {
            return ['KCl', 'H₂O'];
        }
        if (reactants.includes('HCl') && reactants.includes('Ca(OH)₂')) {
            return ['CaCl₂', 'H₂O'];
        }
        if (reactants.includes('H₂SO₄') && reactants.includes('NaOH')) {
            return ['Na₂SO₄', 'H₂O'];
        }
        if (reactants.includes('H₂SO₄') && reactants.includes('KOH')) {
            return ['K₂SO₄', 'H₂O'];
        }
        if (reactants.includes('H₂SO₄') && reactants.includes('Ca(OH)₂')) {
            return ['CaSO₄', 'H₂O'];
        }
        if (reactants.includes('HNO₃') && reactants.includes('NaOH')) {
            return ['NaNO₃', 'H₂O'];
        }
        if (reactants.includes('HNO₃') && reactants.includes('KOH')) {
            return ['KNO₃', 'H₂O'];
        }
        if (reactants.includes('CH₃COOH') && reactants.includes('NaOH')) {
            return ['CH₃COONa', 'H₂O'];
        }

        // Specific acid-base reactions
        if (reactants.includes('SO₂') && reactants.includes('NaOH')) {
            return ['Na₂SO₃', 'H₂O'];
        }
        if (reactants.includes('NO₂') && reactants.includes('H₂O')) {
            return ['HNO₃', 'NO'];
        }
        if (reactants.includes('MgSO₄') && reactants.includes('NaOH')) {
            return ['Mg(OH)₂', 'Na₂SO₄'];
        }
        if (reactants.includes('ZnSO₄') && reactants.includes('NaOH')) {
            return ['Zn(OH)₂', 'Na₂SO₄'];
        }
        if (reactants.includes('FeSO₄') && reactants.includes('NaOH')) {
            return ['Fe(OH)₂', 'Na₂SO₄'];
        }
        if (reactants.includes('AlCl₃') && reactants.includes('NaOH')) {
            return ['Al(OH)₃', 'NaCl'];
        }
        if (reactants.includes('BaCl₂') && reactants.includes('H₂SO₄')) {
            return ['BaSO₄', 'HCl'];
        }
        if (reactants.includes('MgCl₂') && reactants.includes('AgNO₃')) {
            return ['AgCl', 'Mg(NO₃)₂'];
        }

        // Simplified fallback - always produces salt + water
        return ['H₂O']; // Salt would depend on specific acid/base
    }

    // Check if it's a combustion reaction
    isCombustionReaction(reactants) {
        const hydrocarbons = ['CH₄', 'C₂H₆', 'C₃H₈', 'C₄H₁₀', 'C₂H₄', 'C₂H₂', 'C₆H₆', 'C₂H₅OH', 'CH₃OH'];
        const hasHydrocarbon = reactants.some(r => hydrocarbons.includes(r));
        const hasOxygen = reactants.includes('O₂');

        // Must have exactly a hydrocarbon + oxygen (not other compounds)
        return reactants.length === 2 && hasHydrocarbon && hasOxygen;
    }

    // Check if it's metal oxidation
    isMetalOxidation(reactants) {
        const metals = ['Fe', 'Cu', 'Zn', 'Al', 'Mg', 'Ca', 'Na', 'K'];
        const hasOxygen = reactants.includes('O₂');
        const hasMetal = reactants.some(r => metals.includes(r));

        return hasMetal && hasOxygen;
    }

    // Get metal oxidation products
    getMetalOxidationProducts(reactants) {
        const metalOxides = {
            'Fe': 'Fe₂O₃',
            'Cu': 'CuO',
            'Zn': 'ZnO',
            'Al': 'Al₂O₃',
            'Mg': 'MgO',
            'Ca': 'CaO',
            'Na': 'Na₂O',
            'K': 'K₂O'
        };

        const products = [];
        reactants.forEach(r => {
            if (metalOxides[r]) {
                products.push(metalOxides[r]);
            }
        });

        return products;
    }

    // Check if it's double displacement
    isDoubleDisplacement(reactants) {
        if (reactants.length !== 2) return false;

        // Must be two ionic compounds (salts), not metals or simple compounds
        const ionicCompounds = ['AgNO₃', 'NaCl', 'KCl', 'CaCl₂', 'MgCl₂', 'Pb(NO₃)₂', 'KI', 'NaBr', 'KBr', 'BaCl₂', 'Na₂SO₄', 'K₂SO₄', 'Na₂CO₃'];
        const matchingCompounds = reactants.filter(r => ionicCompounds.includes(r));
        console.log('Double displacement check:', reactants, 'matching:', matchingCompounds, 'length:', matchingCompounds.length);

        return matchingCompounds.length === 2;
    }

    // Get double displacement products
    getDoubleDisplacementProducts(reactants) {
        // Common precipitation reactions
        if (reactants.includes('AgNO₃') && reactants.includes('NaCl')) {
            return ['AgCl', 'NaNO₃'];
        }
        if (reactants.includes('Pb(NO₃)₂') && reactants.includes('KI')) {
            return ['PbI₂', 'KNO₃'];
        }
        if (reactants.includes('BaCl₂') && reactants.includes('Na₂SO₄')) {
            return ['BaSO₄', 'NaCl'];
        }
        if (reactants.includes('AgNO₃') && reactants.includes('KBr')) {
            return ['AgBr', 'KNO₃'];
        }
        if (reactants.includes('CaCl₂') && reactants.includes('Na₂CO₃')) {
            return ['CaCO₃', 'NaCl'];
        }
        if (reactants.includes('BaCl₂') && reactants.includes('Na₂SO₄')) {
            return ['BaSO₄', 'NaCl'];
        }
        if (reactants.includes('MgCl₂') && reactants.includes('AgNO₃')) {
            return ['AgCl', 'Mg(NO₃)₂'];
        }
        return [];
    }

    // Show smart suggestions in dropdown
    showSmartSuggestions(type, currentValue) {
        console.log('showSmartSuggestions called with:', type, currentValue);
        this.showingSmartSuggestions = true;

        const compounds = currentValue.split('+').map(c => c.trim()).filter(c => c);
        const suggestions = this.getSmartSuggestions(type, currentValue);
        console.log('Got suggestions:', suggestions);

        const listElement = document.getElementById(`${type}-list`);
        console.log('List element found:', !!listElement);
        if (!listElement) return;

        console.log('Clearing list element and populating with suggestions...');
        listElement.innerHTML = '';

        // If we have smart suggestions, show them first
        if (suggestions.length > 0) {
            console.log('Adding smart suggestions header...');
            // Add header for smart suggestions
            const header = document.createElement('div');
            header.className = 'px-3 py-2 text-xs font-semibold text-gray-500 uppercase tracking-wide suggestion-header';
            header.textContent = '💡 Smart Suggestions';
            listElement.appendChild(header);
            console.log('Header added');

            // Add suggested compounds
            console.log('Adding suggested compounds:', suggestions.slice(0, 8));
            suggestions.slice(0, 8).forEach((compound, index) => {
                console.log(`Processing suggestion ${index + 1}:`, compound);
                const compoundData = this.compounds.find(c => c.formula === compound);
                console.log('Found compound data:', !!compoundData);
                if (compoundData) {
                    try {
                        const option = this.createCompoundOption(compoundData, type);
                        option.classList.add('smart-suggestion');
                        listElement.appendChild(option);
                        console.log(`Added suggestion ${index + 1} to list`);
                    } catch (error) {
                        console.error('Error creating compound option:', error);
                    }
                }
            });
            console.log('All suggestions added');

            // Add separator
            const separator = document.createElement('div');
            separator.className = 'border-t border-gray-200 my-2';
            listElement.appendChild(separator);

            // Add "All Compounds" header
            const allHeader = document.createElement('div');
            allHeader.className = 'px-3 py-2 text-xs font-semibold text-gray-500 uppercase tracking-wide bg-gray-50';
            allHeader.textContent = 'All Compounds';
            listElement.appendChild(allHeader);
        }

        // Group compounds by category
        const categories = {};
        this.compounds.forEach(compound => {
            if (!categories[compound.category]) {
                categories[compound.category] = [];
            }
            categories[compound.category].push(compound);
        });

        // Add all compounds organized by category
        Object.keys(categories).sort().forEach(category => {
            // Add category header
            const categoryHeader = document.createElement('div');
            categoryHeader.className = 'px-3 py-1 text-xs font-medium text-gray-600 bg-gray-50';
            categoryHeader.textContent = category;
            listElement.appendChild(categoryHeader);

            // Add compounds in this category
            categories[category].forEach(compound => {
                // Skip if this compound is already in smart suggestions
                if (suggestions.includes(compound.formula)) return;

                const option = this.createCompoundOption(compound, type);
                listElement.appendChild(option);
            });
        });

        console.log('About to call showDropdown from showSmartSuggestions');
        this.showDropdown(type);
        console.log('showDropdown called, dropdown should be visible now');

        // Double-check the dropdown visibility
        const dropdown = document.getElementById(`${type}-dropdown`);
        if (dropdown) {
            console.log('Dropdown classes after showDropdown:', dropdown.className);
            console.log('Dropdown hidden?', dropdown.classList.contains('hidden'));
        }

        // Reset flag after a short delay
        setTimeout(() => {
            this.showingSmartSuggestions = false;
        }, 100);
    }

    // Create a compound option element
    createCompoundOption(compound, type) {
        const option = document.createElement('div');
        option.className = 'compound-option px-3 py-2 hover:bg-blue-50 cursor-pointer border-b border-gray-100 last:border-b-0';
        option.dataset.formula = compound.formula;
        option.dataset.name = compound.name;

        option.innerHTML = `
            <div class="flex justify-between items-center">
                <span class="font-mono text-lg">${compound.formula}</span>
                <span class="text-sm text-gray-500">${compound.name}</span>
            </div>
        `;

        // Add click handler
        option.addEventListener('click', () => {
            this.selectCompound(type, compound.formula);
        });

        return option;
    }

    // Create a complete equation option element
    createCompleteEquationOption(completeEquation, type, compounds) {
        const option = document.createElement('div');
        option.className = 'compound-option px-3 py-3 hover:bg-blue-50 cursor-pointer border-b border-gray-100 last:border-b-0';
        option.dataset.formula = completeEquation;
        option.dataset.isComplete = 'true';

        // Create a more prominent display for complete equations
        option.innerHTML = `
            <div class="flex flex-col">
                <div class="flex items-center justify-between mb-1">
                    <span class="font-mono text-lg font-semibold text-blue-700">${completeEquation}</span>
                    <span class="text-xs bg-blue-100 text-blue-800 px-2 py-1 rounded-full">Complete</span>
                </div>
                <div class="text-sm text-gray-600">
                    Click to add all products: ${compounds.map(c => this.getCompoundName(c)).join(', ')}
                </div>
            </div>
        `;

        // Add click handler for complete equation
        option.addEventListener('click', () => {
            this.selectCompleteEquation(type, completeEquation);
        });

        return option;
    }

    // Get compound name from formula
    getCompoundName(formula) {
        const compound = this.compounds.find(c => c.formula === formula);
        return compound ? compound.name : formula;
    }

    // Handle selection of complete equation
    selectCompleteEquation(type, completeEquation) {
        const input = document.getElementById(`${type}-input`);
        if (input) {
            input.value = completeEquation;
            input.style.color = ''; // Remove any gray styling
            input.classList.remove('has-suggestion');

            // Close dropdown
            this.closeDropdown(type);

            // Update the other field's placeholder
            if (type === 'products') {
                this.updateReactantPlaceholder();
            } else {
                this.updateProductPlaceholder();
            }
        }
    }

    populateDropdown(type) {
        const listElement = document.getElementById(`${type}-list`);
        if (!listElement) return;

        // Group compounds by category
        const categories = {};
        this.compounds.forEach(compound => {
            if (!categories[compound.category]) {
                categories[compound.category] = [];
            }
            categories[compound.category].push(compound);
        });

        let html = '';
        Object.keys(categories).sort().forEach(category => {
            html += `<div class="px-3 py-2 text-xs font-semibold text-gray-500 uppercase tracking-wide bg-gray-50">${category}</div>`;
            categories[category].forEach(compound => {
                html += `
                    <div class="compound-option px-3 py-2 hover:bg-blue-50 cursor-pointer border-b border-gray-100 last:border-b-0" 
                         data-formula="${compound.formula}" data-name="${compound.name}">
                        <div class="flex justify-between items-center">
                            <span class="font-mono text-lg">${compound.formula}</span>
                            <span class="text-sm text-gray-500">${compound.name}</span>
                        </div>
                    </div>
                `;
            });
        });

        listElement.innerHTML = html;

        // Add click handlers
        listElement.querySelectorAll('.compound-option').forEach(option => {
            option.addEventListener('click', () => {
                this.selectCompound(type, option.dataset.formula);
            });
        });
    }

    filterDropdown(type, query) {
        const listElement = document.getElementById(`${type}-list`);
        if (!listElement) return;

        if (!query.trim()) {
            // Show all compounds organized by category
            this.populateDropdown(type);
            return;
        }

        // Normalize query for better matching (handle subscripts)
        const normalizedQuery = query.toLowerCase()
            .replace(/2/g, '₂').replace(/3/g, '₃').replace(/4/g, '₄')
            .replace(/5/g, '₅').replace(/6/g, '₆').replace(/7/g, '₇')
            .replace(/8/g, '₈').replace(/9/g, '₉').replace(/0/g, '₀');

        // Find matching compounds and score them
        const matches = [];
        this.compounds.forEach(compound => {
            const formula = compound.formula.toLowerCase();
            const name = compound.name.toLowerCase();
            const queryLower = query.toLowerCase();

            let score = 0;

            // Exact formula match gets highest score
            if (formula === queryLower || formula === normalizedQuery.toLowerCase()) {
                score = 1000;
            }
            // Formula starts with query
            else if (formula.startsWith(queryLower) || formula.startsWith(normalizedQuery.toLowerCase())) {
                score = 900;
            }
            // Formula contains query
            else if (formula.includes(queryLower) || formula.includes(normalizedQuery.toLowerCase())) {
                score = 800;
            }
            // Name starts with query
            else if (name.startsWith(queryLower)) {
                score = 700;
            }
            // Name contains query
            else if (name.includes(queryLower)) {
                score = 600;
            }

            // Boost score for common compounds
            const commonCompounds = ['h₂o', 'co₂', 'o₂', 'h₂', 'nacl', 'hcl', 'naoh', 'ch₄', 'nh₃'];
            if (commonCompounds.includes(formula)) {
                score += 100;
            }

            if (score > 0) {
                matches.push({ ...compound, score });
            }
        });

        // Sort by score (highest first)
        matches.sort((a, b) => b.score - a.score);

        // Display filtered results
        let html = '';
        if (matches.length > 0) {
            html += '<div class="px-3 py-2 text-xs font-semibold text-gray-500 uppercase tracking-wide bg-gray-50">Matching Compounds</div>';
            matches.slice(0, 10).forEach(compound => { // Show top 10 matches
                html += `
                    <div class="compound-option px-3 py-2 hover:bg-blue-50 cursor-pointer border-b border-gray-100 last:border-b-0" 
                         data-formula="${compound.formula}" data-name="${compound.name}">
                        <div class="flex justify-between items-center">
                            <span class="font-mono text-lg">${compound.formula}</span>
                            <span class="text-sm text-gray-500">${compound.name}</span>
                        </div>
                    </div>
                `;
            });
        } else {
            html = '<div class="px-3 py-2 text-sm text-gray-500">No compounds found</div>';
        }

        listElement.innerHTML = html;

        // Add click handlers to new options
        listElement.querySelectorAll('.compound-option').forEach(option => {
            option.addEventListener('click', () => {
                this.selectCompound(type, option.dataset.formula);
            });
        });
    }

    selectCompound(type, formula) {
        const input = document.getElementById(`${type}-input`);
        if (!input) return;

        const currentValue = input.value.trim();

        if (currentValue === '') {
            // First compound - just add it
            input.value = formula;
        } else {
            // Check if we're in the middle of typing a compound or adding a new one
            const currentTyping = this.getCurrentCompound(currentValue);

            if (currentTyping && currentTyping.length > 0) {
                // Replace the compound being typed
                const compounds = currentValue.split('+').map(c => c.trim()).filter(c => c.length > 0);
                compounds[compounds.length - 1] = formula;
                input.value = compounds.join(' + ');
            } else {
                // Only add + if the current value doesn't already end with +
                if (currentValue.endsWith(' + ')) {
                    input.value = currentValue + formula;
                } else if (currentValue.endsWith('+')) {
                    input.value = currentValue + ' ' + formula;
                } else {
                    // Add new compound with proper spacing
                    input.value = currentValue + ' + ' + formula;
                }
            }
        }

        input.focus();
        this.closeDropdown(type);

        // Set cursor at the end
        input.setSelectionRange(input.value.length, input.value.length);
    }

    toggleDropdown(type) {
        if (this.activeDropdown === type) {
            this.closeDropdown(type);
        } else {
            this.closeAllDropdowns();
            this.showDropdown(type);
        }
    }

    showDropdown(type) {
        console.log('showDropdown called for:', type);
        const dropdown = document.getElementById(`${type}-dropdown`);
        console.log('Dropdown element found:', !!dropdown);
        if (dropdown) {
            dropdown.classList.remove('hidden');
            this.activeDropdown = type;
            console.log('Dropdown should now be visible');
        }
    }

    closeDropdown(type) {
        console.log('closeDropdown called for:', type);
        const dropdown = document.getElementById(`${type}-dropdown`);
        if (dropdown) {
            dropdown.classList.add('hidden');
            if (this.activeDropdown === type) {
                this.activeDropdown = null;
            }
            console.log('Dropdown closed for:', type);
        }
    }

    closeAllDropdowns() {
        console.log('closeAllDropdowns called');
        this.closeDropdown('reactants');
        this.closeDropdown('products');
    }

    addExamplePlaceholders() {
        // Set initial random examples
        const examples = [
            { reactants: 'CH₄ + O₂', products: 'CO₂ + H₂O' },
            { reactants: 'H₂ + O₂', products: 'H₂O' },
            { reactants: 'Fe + O₂', products: 'Fe₂O₃' },
            { reactants: 'HCl + NaOH', products: 'NaCl + H₂O' },
            { reactants: 'C₂H₆ + O₂', products: 'CO₂ + H₂O' }
        ];

        const reactantsInput = document.getElementById('reactants-input');
        const productsInput = document.getElementById('products-input');

        if (reactantsInput && productsInput) {
            const currentExample = examples[Math.floor(Math.random() * examples.length)];

            // Set initial placeholders
            reactantsInput.placeholder = currentExample.reactants;
            productsInput.placeholder = currentExample.products;

            // Update placeholders dynamically based on input
            this.setupDynamicPlaceholders();
        }
    }

    setupDynamicPlaceholders() {
        const reactantsInput = document.getElementById('reactants-input');
        const productsInput = document.getElementById('products-input');

        if (!reactantsInput || !productsInput) return;

        // Placeholder updates are now handled in the main input event listeners

        // Clear placeholders when focused and empty
        reactantsInput.addEventListener('focus', () => {
            if (!reactantsInput.value.trim()) {
                reactantsInput.placeholder = '';
            }
        });

        productsInput.addEventListener('focus', () => {
            if (!productsInput.value.trim()) {
                productsInput.placeholder = '';
            }
        });

        // Restore default placeholders when fields are empty and blurred
        reactantsInput.addEventListener('blur', () => {
            if (!reactantsInput.value.trim()) {
                reactantsInput.placeholder = 'CH₄ + O₂';
            }
        });

        productsInput.addEventListener('blur', () => {
            if (!productsInput.value.trim()) {
                productsInput.placeholder = 'CO₂ + H₂O';
            }
        });
    }

    // Predict reactants based on products (reverse prediction)
    predictReactants(products) {
        if (products.length === 0) return [];

        // Handle single products
        if (products.length === 1) {
            return this.getSingleProductReactants(products[0]);
        }

        // Handle multiple products with comprehensive reverse patterns
        const reactants = [];

        // Combustion products → Hydrocarbon + O₂
        if (products.includes('CO₂') && products.includes('H₂O')) {
            reactants.push('CH₄', 'O₂');
        }
        // Salt + Water → Acid + Base
        else if (products.includes('H₂O')) {
            if (products.includes('NaCl')) reactants.push('HCl', 'NaOH');
            else if (products.includes('KCl')) reactants.push('HCl', 'KOH');
            else if (products.includes('CaCl₂')) reactants.push('HCl', 'Ca(OH)₂');
            else if (products.includes('MgCl₂')) reactants.push('HCl', 'Mg(OH)₂');
            else if (products.includes('Na₂SO₄')) reactants.push('H₂SO₄', 'NaOH');
            else if (products.includes('K₂SO₄')) reactants.push('H₂SO₄', 'KOH');
            else if (products.includes('CaSO₄')) reactants.push('H₂SO₄', 'Ca(OH)₂');
        }
        // Metal oxide → Metal + O₂
        else if (products.includes('Fe₂O₃')) reactants.push('Fe', 'O₂');
        else if (products.includes('CuO')) reactants.push('Cu', 'O₂');
        else if (products.includes('ZnO')) reactants.push('Zn', 'O₂');
        else if (products.includes('Al₂O₃')) reactants.push('Al', 'O₂');
        else if (products.includes('MgO')) reactants.push('Mg', 'O₂');
        else if (products.includes('CaO')) reactants.push('Ca', 'O₂');

        // Gas + Water → Oxide + Water
        else if (products.includes('CO₂')) {
            if (products.includes('CaO')) reactants.push('CaCO₃');
            else if (products.includes('MgO')) reactants.push('MgCO₃');
            else reactants.push('C', 'O₂');
        }

        // Hydrogen gas production
        else if (products.includes('H₂')) {
            if (products.includes('ZnCl₂')) reactants.push('Zn', 'HCl');
            else if (products.includes('MgCl₂')) reactants.push('Mg', 'HCl');
            else if (products.includes('FeCl₂')) reactants.push('Fe', 'HCl');
            else if (products.includes('NaOH')) reactants.push('Na', 'H₂O');
            else if (products.includes('KOH')) reactants.push('K', 'H₂O');
            else if (products.includes('Ca(OH)₂')) reactants.push('Ca', 'H₂O');
            else reactants.push('H₂', 'O₂');
        }

        // Oxygen gas production
        else if (products.includes('O₂')) {
            if (products.includes('KCl')) reactants.push('KClO₃');
            else if (products.includes('NaCl')) reactants.push('NaClO₃');
            else if (products.includes('H₂O')) reactants.push('H₂O₂');
            else if (products.includes('Ag')) reactants.push('Ag₂O');
            else if (products.includes('Hg')) reactants.push('HgO');
        }

        // Reduction reactions (Metal + CO₂ → Metal oxide + CO)
        else if (products.includes('CO₂')) {
            if (products.includes('Fe')) reactants.push('Fe₂O₃', 'CO');
            else if (products.includes('Cu')) reactants.push('CuO', 'CO');
            else if (products.includes('Pb')) reactants.push('PbO', 'CO');
            else if (products.includes('Mn')) reactants.push('MnO₂', 'C');
        }

        // Halogen displacement
        else if (products.includes('Br₂')) {
            if (products.includes('NaCl')) reactants.push('Cl₂', 'NaBr');
            else if (products.includes('KCl')) reactants.push('Cl₂', 'KBr');
        }
        else if (products.includes('I₂')) {
            if (products.includes('NaCl')) reactants.push('Cl₂', 'NaI');
            else if (products.includes('NaBr')) reactants.push('Br₂', 'NaI');
            else if (products.includes('KCl')) reactants.push('FeCl₃', 'KI');
        }
        else if (products.includes('Cl₂')) {
            if (products.includes('NaF')) reactants.push('F₂', 'NaCl');
            else if (products.includes('KCl')) reactants.push('KMnO₄', 'HCl');
        }

        // Precipitation reactions
        else if (products.includes('AgCl')) {
            if (products.includes('NaNO₃')) reactants.push('AgNO₃', 'NaCl');
        }
        else if (products.includes('PbI₂')) {
            if (products.includes('KNO₃')) reactants.push('Pb(NO₃)₂', 'KI');
        }
        else if (products.includes('AgBr')) {
            if (products.includes('KNO₃')) reactants.push('AgNO₃', 'KBr');
        }

        // Ammonia reactions
        else if (products.includes('NH₄Cl')) {
            reactants.push('NH₃', 'HCl');
        }
        else if (products.includes('(NH₄)₂SO₄')) {
            reactants.push('NH₃', 'H₂SO₄');
        }
        else if (products.includes('NH₃')) {
            if (products.includes('NaCl')) reactants.push('NH₄Cl', 'NaOH');
            else if (products.includes('H₂O')) reactants.push('N₂', 'H₂');
        }

        return reactants;
    }

    // Get reactants for single products
    getSingleProductReactants(product) {
        const singleProductMap = {
            // Simple compounds
            'H₂O': ['H₂', 'O₂'],
            'HCl': ['H₂', 'Cl₂'],
            'HBr': ['H₂', 'Br₂'],
            'HI': ['H₂', 'I₂'],
            'NH₃': ['N₂', 'H₂'],
            'NO': ['N₂', 'O₂'],
            'SO₂': ['S', 'O₂'],
            'CO₂': ['C', 'O₂'],

            // Salts
            'NaCl': ['Na', 'Cl₂'],
            'KCl': ['K', 'Cl₂'],
            'CaCl₂': ['Ca', 'Cl₂'],
            'MgCl₂': ['Mg', 'Cl₂'],
            'FeS': ['Fe', 'S'],
            'CuS': ['Cu', 'S'],
            'ZnS': ['Zn', 'S'],

            // Metal oxides
            'Fe₂O₃': ['Fe', 'O₂'],
            'CuO': ['Cu', 'O₂'],
            'ZnO': ['Zn', 'O₂'],
            'Al₂O₃': ['Al', 'O₂'],
            'MgO': ['Mg', 'O₂'],
            'CaO': ['Ca', 'O₂'],
            'Na₂O': ['Na', 'O₂'],
            'K₂O': ['K', 'O₂'],

            // Hydroxides
            'NaOH': ['Na', 'H₂O'],
            'KOH': ['K', 'H₂O'],
            'Ca(OH)₂': ['Ca', 'H₂O'],
            'Mg(OH)₂': ['Mg', 'H₂O'],

            // Acids from oxides
            'H₂SO₄': ['SO₃', 'H₂O'],
            'H₂CO₃': ['CO₂', 'H₂O'],
            'H₃PO₄': ['P₂O₅', 'H₂O'],
        };

        return singleProductMap[product] || [];
    }

    balanceEquation() {
        console.log('balanceEquation called');

        const reactantsInput = document.getElementById('reactants-input');
        const productsInput = document.getElementById('products-input');

        if (!reactantsInput || !productsInput) {
            console.error('Input elements not found');
            this.showError('Input fields not found');
            return;
        }

        let reactants = reactantsInput.value.trim();
        let products = productsInput.value.trim();

        // Remove trailing plus signs and clean up
        reactants = reactants.replace(/\+\s*$/, '').trim();
        products = products.replace(/\+\s*$/, '').trim();

        // Update the input fields to show cleaned values
        reactantsInput.value = reactants;
        productsInput.value = products;

        console.log('Reactants:', reactants);
        console.log('Products:', products);

        if (!reactants || !products) {
            this.showError('Please enter both reactants and products');
            return;
        }

        // Check for empty compounds (double plus signs, etc.)
        if (reactants.includes('++') || products.includes('++') ||
            reactants.startsWith('+') || products.startsWith('+')) {
            this.showError('Invalid format: Please check for extra plus signs or empty compounds');
            return;
        }

        // Combine into a single equation format
        const equation = `${reactants} -> ${products}`;
        console.log('Combined equation:', equation);

        this.hideError();
        this.hideResults();

        const result = this.balancer.balance(equation);
        console.log('Balance result:', result);

        if (result.success) {
            // Check if equation is already balanced
            if (this.isAlreadyBalanced(result)) {
                this.showAlreadyBalanced(result, equation);
            } else {
                this.showResults(result, equation);
            }
        } else {
            this.showSmartError(result.error, equation);
        }
    }

    // Check if the equation is already balanced (all coefficients are 1)
    isAlreadyBalanced(result) {
        return result.coefficients.every(coeff => coeff === 1);
    }

    // Show message for already balanced equations
    showAlreadyBalanced(result, originalEquation) {
        const resultsSection = document.getElementById('results-section');
        const balancedEquationDiv = document.getElementById('balanced-equation');
        const stepsContainer = document.getElementById('steps-container');

        if (!resultsSection || !balancedEquationDiv || !stepsContainer) {
            console.error('Results elements not found');
            return;
        }

        // Show the equation
        balancedEquationDiv.textContent = result.balancedEquation;

        // Show explanation for already balanced equation
        stepsContainer.innerHTML = `
            <div class="bg-green-50 border border-green-200 rounded-lg p-4">
                <div class="flex items-center mb-2">
                    <i class="fas fa-check-circle text-green-600 mr-2"></i>
                    <h3 class="text-lg font-semibold text-green-800">Equation Already Balanced!</h3>
                </div>
                <p class="text-green-700 mb-3">
                    This equation is already balanced as written. All elements have the same number of atoms on both sides.
                </p>
                <div class="bg-white p-3 rounded border">
                    <h4 class="font-medium text-gray-800 mb-2">Element Verification:</h4>
                    ${this.generateElementVerification(result)}
                </div>
            </div>
        `;

        resultsSection.classList.remove('hidden');
        resultsSection.scrollIntoView({ behavior: 'smooth' });
    }

    // Generate element verification for already balanced equations
    generateElementVerification(result) {
        let content = '';
        result.elements.forEach(element => {
            let reactantCount = 0;
            let productCount = 0;

            result.compounds.forEach((compound, i) => {
                const count = (compound.elements[element] || 0) * result.coefficients[i];
                if (compound.isReactant) {
                    reactantCount += count;
                } else {
                    productCount += count;
                }
            });

            content += `<p class="text-sm"><strong>${element}:</strong> ${reactantCount} (reactants) = ${productCount} (products) ✓</p>`;
        });
        return content;
    }

    // Show smarter error messages
    showSmartError(errorMessage, equation) {
        let smartMessage = errorMessage;
        let suggestions = [];

        // Analyze the error and provide helpful suggestions
        if (errorMessage.includes('Could not find a valid solution')) {
            smartMessage = 'This equation cannot be balanced';
            suggestions = [
                'Check if all compounds are written correctly',
                'Verify that this is a valid chemical reaction',
                'Some reactions may require additional reactants or products',
                'Consider if this might be a nuclear reaction (not supported)'
            ];
        } else if (errorMessage.includes('Empty compound')) {
            smartMessage = 'Empty compound detected';
            suggestions = [
                'Check for extra plus signs (+ +)',
                'Make sure all compounds are properly entered',
                'Remove any trailing plus signs'
            ];
        } else if (errorMessage.includes('Could not parse compound')) {
            smartMessage = 'Invalid compound formula detected';
            suggestions = [
                'Check chemical formulas for typos',
                'Use proper capitalization (H2O, not h2o)',
                'Make sure all elements are valid',
                'Use numbers for subscripts (H2O, not H₂O if typing manually)'
            ];
        } else if (errorMessage.includes('Invalid equation format')) {
            smartMessage = 'Equation format error';
            suggestions = [
                'Make sure both reactants and products are filled in',
                'Check for proper compound separation with + signs',
                'Verify all chemical formulas are complete'
            ];
        }

        // Display the smart error
        const errorSection = document.getElementById('error-section');
        const errorMessage_elem = document.getElementById('error-message');

        if (errorSection && errorMessage_elem) {
            let content = `<strong>${smartMessage}</strong>`;

            if (suggestions.length > 0) {
                content += '<br><br><strong>Suggestions:</strong><ul class="mt-2 ml-4">';
                suggestions.forEach(suggestion => {
                    content += `<li class="text-sm">• ${suggestion}</li>`;
                });
                content += '</ul>';
            }

            content += `<br><br><em class="text-sm">Original equation: ${equation}</em>`;

            errorMessage_elem.innerHTML = content;
            errorSection.classList.remove('hidden');
            errorSection.scrollIntoView({ behavior: 'smooth' });
        }
    }

    showResults(result, originalEquation) {
        console.log('Showing results:', result);

        const resultsSection = document.getElementById('results-section');
        const balancedEquationDiv = document.getElementById('balanced-equation');
        const stepsContainer = document.getElementById('steps-container');

        if (!resultsSection || !balancedEquationDiv || !stepsContainer) {
            console.error('Results elements not found');
            return;
        }

        // Show balanced equation
        balancedEquationDiv.textContent = result.balancedEquation;

        // Generate step-by-step explanation
        this.generateSteps(stepsContainer, result, originalEquation);

        resultsSection.classList.remove('hidden');
        resultsSection.scrollIntoView({ behavior: 'smooth' });
    }

    generateSteps(container, result, originalEquation) {
        container.innerHTML = '';

        const steps = [
            {
                title: '🔍 Step 1: Analyze the Equation & Count Atoms',
                icon: 'fas fa-search',
                color: 'blue',
                content: this.generateAnalyzeStep(originalEquation, result)
            },
            {
                title: '⚖️ Step 2: Balance Each Element & Watch the Equation Change',
                icon: 'fas fa-balance-scale',
                color: 'indigo',
                content: this.generateDynamicBalancingStep(result, originalEquation)
            },
            {
                title: '✅ Step 3: Verify the Final Balance',
                icon: 'fas fa-check-circle',
                color: 'green',
                content: this.generateVerificationStep(result)
            }
        ];

        steps.forEach((step, index) => {
            const stepDiv = document.createElement('div');
            stepDiv.className = `mb-6 border-2 border-${step.color}-200 rounded-xl shadow-lg overflow-hidden`;

            stepDiv.innerHTML = `
                <button class="w-full text-left p-5 bg-gradient-to-r from-${step.color}-50 to-${step.color}-100 hover:from-${step.color}-100 hover:to-${step.color}-150 transition-all duration-300 flex items-center justify-between group" 
                        onclick="this.nextElementSibling.classList.toggle('hidden'); this.querySelector('.chevron').classList.toggle('rotate-180')">
                    <div class="flex items-center">
                        <div class="w-10 h-10 bg-${step.color}-500 text-white rounded-full flex items-center justify-center mr-4 group-hover:scale-110 transition-transform">
                            <i class="${step.icon} text-sm"></i>
                        </div>
                        <span class="font-bold text-lg text-${step.color}-800">${step.title}</span>
                    </div>
                    <i class="fas fa-chevron-down text-${step.color}-600 chevron transition-transform duration-300"></i>
                </button>
                <div class="p-6 bg-white border-t-2 border-${step.color}-100 ${index === 0 ? '' : 'hidden'}" style="position: relative; overflow: hidden; contain: layout;">
                    <div style="position: relative; z-index: 1;">
                        ${step.content}
                    </div>
                </div>
            `;

            container.appendChild(stepDiv);
        });
    }

    generateIdentifyStep(originalEquation, result) {
        const reactants = result.compounds.filter(c => c.isReactant);
        const products = result.compounds.filter(c => !c.isReactant);

        return `
            <div class="space-y-4">
                <div class="bg-gradient-to-r from-blue-50 to-blue-100 p-4 rounded-lg border-l-4 border-blue-500">
                    <h4 class="font-bold text-blue-800 mb-2 flex items-center">
                        <i class="fas fa-flask mr-2"></i>Our Unbalanced Equation
                    </h4>
                    <div class="text-xl font-mono bg-white p-3 rounded border text-center">
                        ${originalEquation}
                    </div>
                    <p class="text-blue-700 text-sm mt-2 text-center">
                        <i class="fas fa-exclamation-triangle mr-1"></i>
                        This equation is NOT balanced yet - we need to add coefficients!
                    </p>
                </div>
                
                <div class="grid md:grid-cols-2 gap-4">
                    <div class="bg-red-50 p-4 rounded-lg border-l-4 border-red-400">
                        <h4 class="font-bold text-red-700 mb-3 flex items-center">
                            <i class="fas fa-play mr-2"></i>What We Start With
                        </h4>
                        <div class="space-y-2">
                            ${reactants.map(c => `
                                <div class="bg-white p-2 rounded border flex items-center">
                                    <span class="font-mono text-lg font-bold text-red-600">${c.formula}</span>
                                </div>
                            `).join('')}
                        </div>
                    </div>
                    
                    <div class="bg-green-50 p-4 rounded-lg border-l-4 border-green-400">
                        <h4 class="font-bold text-green-700 mb-3 flex items-center">
                            <i class="fas fa-flag-checkered mr-2"></i>What We End Up With
                        </h4>
                        <div class="space-y-2">
                            ${products.map(c => `
                                <div class="bg-white p-2 rounded border flex items-center">
                                    <span class="font-mono text-lg font-bold text-green-600">${c.formula}</span>
                                </div>
                            `).join('')}
                        </div>
                    </div>
                </div>
                
                <div class="bg-blue-50 p-4 rounded-lg border border-blue-200">
                    <div class="flex items-start">
                        <i class="fas fa-lightbulb text-blue-600 mt-1 mr-3"></i>
                        <div>
                            <h5 class="font-semibold text-blue-800 mb-1">The Goal</h5>
                            <p class="text-blue-700 text-sm">
                                We need to find the right numbers (coefficients) to put in front of each compound so that 
                                we have the same number of each type of atom on both sides. Think of it like a recipe - 
                                we need the right amounts of each ingredient!
                            </p>
                        </div>
                    </div>
                </div>
            </div>
        `;
    }

    generateElementsStep(result) {
        const elementColors = {
            'H': 'bg-red-100 text-red-800 border-red-300',
            'O': 'bg-blue-100 text-blue-800 border-blue-300',
            'C': 'bg-gray-100 text-gray-800 border-gray-300',
            'N': 'bg-green-100 text-green-800 border-green-300',
            'S': 'bg-yellow-100 text-yellow-800 border-yellow-300',
            'Cl': 'bg-purple-100 text-purple-800 border-purple-300',
            'Na': 'bg-orange-100 text-orange-800 border-orange-300',
            'Ca': 'bg-pink-100 text-pink-800 border-pink-300',
            'Fe': 'bg-indigo-100 text-indigo-800 border-indigo-300'
        };

        return `
            <div class="space-y-4">
                <div class="bg-gradient-to-r from-purple-50 to-purple-100 p-4 rounded-lg border-l-4 border-purple-500">
                    <h4 class="font-bold text-purple-800 mb-3 flex items-center">
                        <i class="fas fa-atom mr-2"></i>Elements Discovered
                    </h4>
                    <div class="flex flex-wrap gap-2">
                        ${result.elements.map(element => {
            const colorClass = elementColors[element] || 'bg-gray-100 text-gray-800 border-gray-300';
            return `
                                <span class="px-3 py-1 rounded-full border-2 font-bold ${colorClass}">
                                    ${element}
                                </span>
                            `;
        }).join('')}
                    </div>
                </div>
                
                <div class="bg-white border-2 border-purple-200 rounded-lg p-4">
                    <h4 class="font-bold text-purple-800 mb-4 flex items-center">
                        <i class="fas fa-microscope mr-2"></i>Element Analysis by Compound
                    </h4>
                    <div class="grid gap-3">
                        ${result.compounds.map(compound => {
            const isReactant = compound.isReactant;
            const bgColor = isReactant ? 'bg-red-50 border-red-200' : 'bg-green-50 border-green-200';
            const textColor = isReactant ? 'text-red-700' : 'text-green-700';
            const icon = isReactant ? 'fas fa-arrow-right' : 'fas fa-bullseye';

            return `
                                <div class="${bgColor} border-2 rounded-lg p-3">
                                    <div class="flex items-center mb-2">
                                        <i class="${icon} ${textColor} mr-2"></i>
                                        <span class="font-mono text-lg font-bold ${textColor}">${compound.formula}</span>
                                        <span class="ml-2 text-sm ${textColor}">(${isReactant ? 'Reactant' : 'Product'})</span>
                                    </div>
                                    <div class="flex flex-wrap gap-2">
                                        ${Object.entries(compound.elements).map(([element, count]) => {
                const colorClass = elementColors[element] || 'bg-gray-100 text-gray-800 border-gray-300';
                return `
                                                <div class="flex items-center bg-white rounded border px-2 py-1">
                                                    <span class="w-6 h-6 rounded-full border-2 ${colorClass} flex items-center justify-center text-xs font-bold mr-2">
                                                        ${element}
                                                    </span>
                                                    <span class="text-sm font-semibold">× ${count}</span>
                                                </div>
                                            `;
            }).join('')}
                                    </div>
                                </div>
                            `;
        }).join('')}
                    </div>
                </div>
                
                <div class="bg-purple-50 p-4 rounded-lg border border-purple-200">
                    <div class="flex items-start">
                        <i class="fas fa-lightbulb text-purple-600 mt-1 mr-3"></i>
                        <div>
                            <h5 class="font-semibold text-purple-800 mb-1">Why do we do this?</h5>
                            <p class="text-purple-700 text-sm">
                                By counting each type of atom in every compound, we can set up equations to ensure the same number 
                                of each element appears on both sides of the reaction. This follows the Law of Conservation of Mass.
                            </p>
                        </div>
                    </div>
                </div>
            </div>
        `;
    }

    generateAnalyzeStep(originalEquation, result) {
        const elementColors = {
            'H': 'bg-red-100 text-red-800 border-red-300',
            'O': 'bg-blue-100 text-blue-800 border-blue-300',
            'C': 'bg-gray-100 text-gray-800 border-gray-300',
            'N': 'bg-green-100 text-green-800 border-green-300',
            'S': 'bg-yellow-100 text-yellow-800 border-yellow-300',
            'Cl': 'bg-purple-100 text-purple-800 border-purple-300',
            'Na': 'bg-orange-100 text-orange-800 border-orange-300',
            'Ca': 'bg-pink-100 text-pink-800 border-pink-300',
            'Fe': 'bg-indigo-100 text-indigo-800 border-indigo-300'
        };

        const reactants = result.compounds.filter(c => c.isReactant);
        const products = result.compounds.filter(c => !c.isReactant);

        // Count atoms before balancing (assuming coefficient of 1 for all)
        const atomCounts = {};
        result.elements.forEach(element => {
            let reactantCount = 0;
            let productCount = 0;

            result.compounds.forEach(compound => {
                const count = compound.elements[element] || 0;
                if (compound.isReactant) {
                    reactantCount += count;
                } else {
                    productCount += count;
                }
            });

            atomCounts[element] = { reactants: reactantCount, products: productCount };
        });

        return `
            <div class="space-y-4">
                <div class="bg-gradient-to-r from-blue-50 to-blue-100 p-4 rounded-lg border-l-4 border-blue-500">
                    <h4 class="font-bold text-blue-800 mb-2 flex items-center">
                        <i class="fas fa-flask mr-2"></i>Our Unbalanced Equation
                    </h4>
                    <div class="text-xl font-mono bg-white p-3 rounded border text-center">
                        ${originalEquation}
                    </div>
                    <p class="text-blue-700 text-sm mt-2 text-center">
                        <i class="fas fa-exclamation-triangle mr-1"></i>
                        This equation is NOT balanced yet - let's see why!
                    </p>
                </div>
                
                <div class="grid md:grid-cols-2 gap-4">
                    <div class="bg-red-50 p-4 rounded-lg border-l-4 border-red-400">
                        <h4 class="font-bold text-red-700 mb-3 flex items-center">
                            <i class="fas fa-play mr-2"></i>Reactants (Left Side)
                        </h4>
                        <div class="space-y-2">
                            ${reactants.map(c => `
                                <div class="bg-white p-2 rounded border">
                                    <span class="font-mono text-lg font-bold text-red-600">${c.formula}</span>
                                    <div class="text-xs text-red-600 mt-1">
                                        ${Object.entries(c.elements).map(([el, count]) => `${el}: ${count}`).join(', ')}
                                    </div>
                                </div>
                            `).join('')}
                        </div>
                    </div>
                    
                    <div class="bg-green-50 p-4 rounded-lg border-l-4 border-green-400">
                        <h4 class="font-bold text-green-700 mb-3 flex items-center">
                            <i class="fas fa-flag-checkered mr-2"></i>Products (Right Side)
                        </h4>
                        <div class="space-y-2">
                            ${products.map(c => `
                                <div class="bg-white p-2 rounded border">
                                    <span class="font-mono text-lg font-bold text-green-600">${c.formula}</span>
                                    <div class="text-xs text-green-600 mt-1">
                                        ${Object.entries(c.elements).map(([el, count]) => `${el}: ${count}`).join(', ')}
                                    </div>
                                </div>
                            `).join('')}
                        </div>
                    </div>
                </div>

                <div class="bg-white border-2 border-blue-200 rounded-lg overflow-hidden">
                    <div class="bg-blue-100 p-3 border-b border-blue-200">
                        <h4 class="font-bold text-blue-800 flex items-center">
                            <i class="fas fa-table mr-2"></i>Atom Count T-Table
                        </h4>
                        <p class="text-blue-700 text-sm mt-1">Let's count atoms on each side to see the imbalance</p>
                    </div>
                    <div class="overflow-x-auto">
                        <table class="w-full">
                            <thead class="bg-gray-50">
                                <tr>
                                    <th class="text-left p-4 font-bold text-gray-700 border-r">Element</th>
                                    <th class="text-center p-4 font-bold text-blue-700 bg-blue-50 border-r">Left Side<br><span class="text-xs font-normal">(Reactants)</span></th>
                                    <th class="text-center p-4 font-bold text-purple-700 bg-purple-50 border-r">Right Side<br><span class="text-xs font-normal">(Products)</span></th>
                                    <th class="text-center p-4 font-bold text-gray-700">Balanced?</th>
                                </tr>
                            </thead>
                            <tbody>
                                ${result.elements.map(element => {
            const colorClass = elementColors[element] || 'bg-gray-100 text-gray-800 border-gray-300';
            const counts = atomCounts[element];
            const isBalanced = counts.reactants === counts.products;

            return `
                                        <tr class="border-t hover:bg-gray-50">
                                            <td class="p-4 border-r">
                                                <div class="flex items-center">
                                                    <span class="w-8 h-8 rounded-full border-2 ${colorClass} flex items-center justify-center text-sm font-bold mr-3">
                                                        ${element}
                                                    </span>
                                                    <span class="font-bold text-gray-800">${element}</span>
                                                </div>
                                            </td>
                                            <td class="text-center p-4 bg-blue-50 border-r">
                                                <div class="text-2xl font-bold text-blue-600">${counts.reactants}</div>
                                            </td>
                                            <td class="text-center p-4 bg-purple-50 border-r">
                                                <div class="text-2xl font-bold text-purple-600">${counts.products}</div>
                                            </td>
                                            <td class="text-center p-4">
                                                <div class="text-2xl ${isBalanced ? 'text-green-500' : 'text-red-500'}">
                                                    <i class="fas fa-${isBalanced ? 'check-circle' : 'times-circle'}"></i>
                                                </div>
                                                <div class="text-xs ${isBalanced ? 'text-green-600' : 'text-red-600'} mt-1">
                                                    ${isBalanced ? 'Balanced' : 'Unbalanced'}
                                                </div>
                                            </td>
                                        </tr>
                                    `;
        }).join('')}
                            </tbody>
                        </table>
                    </div>
                </div>
                
                <div class="bg-blue-50 p-4 rounded-lg border border-blue-200">
                    <div class="flex items-start">
                        <i class="fas fa-lightbulb text-blue-600 mt-1 mr-3"></i>
                        <div>
                            <h5 class="font-semibold text-blue-800 mb-1">The Problem is Clear!</h5>
                            <p class="text-blue-700 text-sm">
                                The T-table shows us exactly which elements are unbalanced. We need to add coefficients 
                                (numbers in front of compounds) to make the left and right sides equal for every element.
                            </p>
                        </div>
                    </div>
                </div>
            </div>
        `;
    }

    generateCountStep(result) {
        const elementColors = {
            'H': 'bg-red-100 text-red-800 border-red-300',
            'O': 'bg-blue-100 text-blue-800 border-blue-300',
            'C': 'bg-gray-100 text-gray-800 border-gray-300',
            'N': 'bg-green-100 text-green-800 border-green-300',
            'S': 'bg-yellow-100 text-yellow-800 border-yellow-300',
            'Cl': 'bg-purple-100 text-purple-800 border-purple-300',
            'Na': 'bg-orange-100 text-orange-800 border-orange-300',
            'Ca': 'bg-pink-100 text-pink-800 border-pink-300',
            'Fe': 'bg-indigo-100 text-indigo-800 border-indigo-300'
        };

        // Count atoms before balancing (assuming coefficient of 1 for all)
        const atomCounts = {};
        result.elements.forEach(element => {
            let reactantCount = 0;
            let productCount = 0;

            result.compounds.forEach(compound => {
                const count = compound.elements[element] || 0;
                if (compound.isReactant) {
                    reactantCount += count;
                } else {
                    productCount += count;
                }
            });

            atomCounts[element] = { reactants: reactantCount, products: productCount };
        });

        return `
            <div class="space-y-4">
                <div class="bg-gradient-to-r from-purple-50 to-purple-100 p-4 rounded-lg border-l-4 border-purple-500">
                    <h4 class="font-bold text-purple-800 mb-2 flex items-center">
                        <i class="fas fa-calculator mr-2"></i>Let's Count the Atoms
                    </h4>
                    <p class="text-purple-700 text-sm">
                        Before we balance, let's see what we have. We'll count each type of atom on both sides.
                    </p>
                </div>
                
                <div class="bg-white border-2 border-purple-200 rounded-lg p-4">
                    <h4 class="font-bold text-purple-800 mb-4 flex items-center">
                        <i class="fas fa-microscope mr-2"></i>Atom Inventory (Before Balancing)
                    </h4>
                    <div class="space-y-3">
                        ${result.elements.map(element => {
            const colorClass = elementColors[element] || 'bg-gray-100 text-gray-800 border-gray-300';
            const counts = atomCounts[element];
            const isBalanced = counts.reactants === counts.products;

            return `
                                <div class="bg-gray-50 border-2 border-gray-200 rounded-lg p-4">
                                    <div class="flex items-center justify-between">
                                        <div class="flex items-center">
                                            <span class="w-8 h-8 rounded-full border-2 ${colorClass} flex items-center justify-center text-sm font-bold mr-3">
                                                ${element}
                                            </span>
                                            <span class="font-bold text-gray-800">${element} atoms:</span>
                                        </div>
                                        <div class="flex items-center space-x-4">
                                            <div class="text-center">
                                                <div class="bg-red-100 text-red-800 px-3 py-1 rounded-full font-bold">
                                                    ${counts.reactants}
                                                </div>
                                                <div class="text-xs text-red-600 mt-1">Left Side</div>
                                            </div>
                                            <div class="text-2xl font-bold ${isBalanced ? 'text-green-500' : 'text-red-500'}">
                                                ${isBalanced ? '=' : '≠'}
                                            </div>
                                            <div class="text-center">
                                                <div class="bg-green-100 text-green-800 px-3 py-1 rounded-full font-bold">
                                                    ${counts.products}
                                                </div>
                                                <div class="text-xs text-green-600 mt-1">Right Side</div>
                                            </div>
                                            <div class="text-lg ${isBalanced ? 'text-green-500' : 'text-red-500'}">
                                                <i class="fas fa-${isBalanced ? 'check' : 'times'}"></i>
                                            </div>
                                        </div>
                                    </div>
                                </div>
                            `;
        }).join('')}
                    </div>
                </div>
                
                <div class="bg-purple-50 p-4 rounded-lg border border-purple-200">
                    <div class="flex items-start">
                        <i class="fas fa-exclamation-triangle text-purple-600 mt-1 mr-3"></i>
                        <div>
                            <h5 class="font-semibold text-purple-800 mb-1">The Problem</h5>
                            <p class="text-purple-700 text-sm">
                                As you can see, most elements don't have equal numbers on both sides. That's why we need to balance! 
                                We'll add coefficients (numbers in front of compounds) to make both sides equal.
                            </p>
                        </div>
                    </div>
                </div>
            </div>
        `;
    }

    generateBalancingStep(result, originalEquation) {
        // Determine balancing strategy based on elements present
        const hasMetals = result.elements.some(el => ['Na', 'K', 'Ca', 'Mg', 'Fe', 'Cu', 'Zn', 'Al', 'Ag', 'Pb', 'Sn', 'Cr', 'Mn'].includes(el));
        const hasHydrogen = result.elements.includes('H');
        const hasOxygen = result.elements.includes('O');
        const hasCarbon = result.elements.includes('C');

        // Create balancing order based on best practices
        let balancingOrder = [];
        let strategy = "";

        if (hasMetals) {
            strategy = "We'll balance metals first (they appear in fewer compounds), then nonmetals, hydrogen second-to-last, and oxygen last.";
            balancingOrder = result.elements.filter(el => ['Na', 'K', 'Ca', 'Mg', 'Fe', 'Cu', 'Zn', 'Al', 'Ag', 'Pb', 'Sn', 'Cr', 'Mn'].includes(el));
            balancingOrder = balancingOrder.concat(result.elements.filter(el => !['Na', 'K', 'Ca', 'Mg', 'Fe', 'Cu', 'Zn', 'Al', 'Ag', 'Pb', 'Sn', 'Cr', 'Mn', 'H', 'O'].includes(el)));
            if (hasHydrogen) balancingOrder.push('H');
            if (hasOxygen) balancingOrder.push('O');
        } else if (hasCarbon) {
            strategy = "For organic reactions: carbon first (molecular backbone), then hydrogen, then oxygen last.";
            balancingOrder = ['C'];
            if (hasHydrogen) balancingOrder.push('H');
            if (hasOxygen) balancingOrder.push('O');
            balancingOrder = balancingOrder.concat(result.elements.filter(el => !['C', 'H', 'O'].includes(el)));
        } else {
            strategy = "We'll balance elements that appear in fewer compounds first, saving polyatomic elements for last.";
            balancingOrder = [...result.elements];
        }

        return `
            <div class="space-y-4">
                <div class="bg-gradient-to-r from-indigo-50 to-indigo-100 p-4 rounded-lg border-l-4 border-indigo-500">
                    <h4 class="font-bold text-indigo-800 mb-2 flex items-center">
                        <i class="fas fa-route mr-2"></i>Why Order Matters in Balancing
                    </h4>
                    <p class="text-indigo-700 text-sm mb-3">${strategy}</p>
                    <div class="bg-white p-3 rounded border">
                        <h5 class="font-semibold text-indigo-800 text-sm mb-2">The Logic:</h5>
                        <ul class="text-xs text-indigo-700 space-y-1">
                            <li><strong>Metals first:</strong> Usually appear in only 1-2 compounds, so easier to balance</li>
                            <li><strong>Nonmetals next:</strong> More complex but still manageable</li>
                            <li><strong>Hydrogen second-to-last:</strong> Often appears in many compounds (acids, bases, water)</li>
                            <li><strong>Oxygen last:</strong> Appears in the most compounds, hardest to balance early</li>
                        </ul>
                    </div>
                </div>
                
                <div class="bg-white border-2 border-indigo-200 rounded-lg p-4">
                    <h4 class="font-bold text-indigo-800 mb-4 flex items-center">
                        <i class="fas fa-step-forward mr-2"></i>How We Balance Each Element
                    </h4>
                    <div class="space-y-4">
                        ${balancingOrder.map((element, index) => {
            const elementColors = {
                'H': 'bg-pink-100 text-pink-800 border-pink-300',
                'O': 'bg-blue-100 text-blue-800 border-blue-300',
                'C': 'bg-gray-100 text-gray-800 border-gray-300',
                'N': 'bg-teal-100 text-teal-800 border-teal-300',
                'S': 'bg-yellow-100 text-yellow-800 border-yellow-300',
                'Cl': 'bg-purple-100 text-purple-800 border-purple-300',
                'Na': 'bg-orange-100 text-orange-800 border-orange-300',
                'Ca': 'bg-lime-100 text-lime-800 border-lime-300',
                'Fe': 'bg-amber-100 text-amber-800 border-amber-300'
            };
            const colorClass = elementColors[element] || 'bg-slate-100 text-slate-800 border-slate-300';

            // Calculate current atom counts for this element
            let leftSideCount = 0;
            let rightSideCount = 0;
            let leftCompounds = [];
            let rightCompounds = [];

            result.compounds.forEach((compound, i) => {
                const count = (compound.elements[element] || 0) * result.coefficients[i];
                if (compound.isReactant) {
                    leftSideCount += count;
                    if (compound.elements[element]) {
                        leftCompounds.push(`${result.coefficients[i] === 1 ? '' : result.coefficients[i]}${compound.formula}`);
                    }
                } else {
                    rightSideCount += count;
                    if (compound.elements[element]) {
                        rightCompounds.push(`${result.coefficients[i] === 1 ? '' : result.coefficients[i]}${compound.formula}`);
                    }
                }
            });

            let balancingExplanation = "";
            if (['Na', 'K', 'Ca', 'Mg', 'Fe', 'Cu', 'Zn', 'Al', 'Ag', 'Pb', 'Sn', 'Cr', 'Mn'].includes(element)) {
                balancingExplanation = `Metal atoms are usually the easiest to balance because they appear in fewer compounds. We look for the compound with the most ${element} atoms and work from there.`;
            } else if (element === 'C') {
                balancingExplanation = "Carbon forms the backbone of organic molecules. Balance it first to establish the molecular ratios.";
            } else if (element === 'H') {
                balancingExplanation = "Hydrogen appears in many compounds (water, acids, bases). We balance it after other elements are set.";
            } else if (element === 'O') {
                balancingExplanation = "Oxygen appears in the most compounds and is saved for last. By this point, the coefficients usually balance oxygen automatically.";
            } else {
                balancingExplanation = `This nonmetal is balanced after metals but before hydrogen and oxygen.`;
            }

            return `
                                <div class="border-2 border-gray-200 rounded-lg p-4">
                                    <div class="flex items-center mb-3">
                                        <div class="w-8 h-8 bg-indigo-500 text-white rounded-full flex items-center justify-center font-bold mr-3">
                                            ${index + 1}
                                        </div>
                                        <span class="w-10 h-10 rounded-full border-2 ${colorClass} flex items-center justify-center text-sm font-bold mr-3">
                                            ${element}
                                        </span>
                                        <div>
                                            <div class="font-bold text-gray-800">Balancing ${element}</div>
                                            <div class="text-sm text-gray-600">${balancingExplanation}</div>
                                        </div>
                                    </div>
                                    
                                    <div class="bg-gray-50 p-3 rounded-lg">
                                        <div class="grid grid-cols-3 gap-4 items-center">
                                            <div class="text-center">
                                                <div class="bg-blue-100 text-blue-800 px-3 py-2 rounded-lg font-bold text-lg">
                                                    ${leftSideCount}
                                                </div>
                                                <div class="text-xs text-blue-600 mt-1">Left Side</div>
                                                <div class="text-xs text-gray-600 mt-1">${leftCompounds.join(' + ')}</div>
                                            </div>
                                            
                                            <div class="text-center">
                                                <div class="text-2xl font-bold ${leftSideCount === rightSideCount ? 'text-green-600' : 'text-orange-600'}">
                                                    ${leftSideCount === rightSideCount ? '=' : '→'}
                                                </div>
                                                <div class="text-xs ${leftSideCount === rightSideCount ? 'text-green-600' : 'text-orange-600'} mt-1">
                                                    ${leftSideCount === rightSideCount ? 'Balanced!' : 'Balancing...'}
                                                </div>
                                            </div>
                                            
                                            <div class="text-center">
                                                <div class="bg-purple-100 text-purple-800 px-3 py-2 rounded-lg font-bold text-lg">
                                                    ${rightSideCount}
                                                </div>
                                                <div class="text-xs text-purple-600 mt-1">Right Side</div>
                                                <div class="text-xs text-gray-600 mt-1">${rightCompounds.join(' + ')}</div>
                                            </div>
                                        </div>
                                    </div>
                                </div>
                            `;
        }).join('')}
                    </div>
                </div>
                
                <div class="bg-indigo-50 p-4 rounded-lg border border-indigo-200">
                    <div class="flex items-start">
                        <i class="fas fa-brain text-indigo-600 mt-1 mr-3"></i>
                        <div>
                            <h5 class="font-semibold text-indigo-800 mb-1">The Strategy in Action</h5>
                            <p class="text-indigo-700 text-sm">
                                Notice how we work systematically through each element. By following this order, 
                                we avoid the frustration of constantly changing coefficients. Each element builds 
                                on the previous ones until the entire equation is balanced!
                            </p>
                        </div>
                    </div>
                </div>
            </div>
        `;
    }

    generateDynamicBalancingStep(result, originalEquation) {
        // Determine balancing strategy
        const hasMetals = result.elements.some(el => ['Na', 'K', 'Ca', 'Mg', 'Fe', 'Cu', 'Zn', 'Al', 'Ag', 'Pb', 'Sn', 'Cr', 'Mn'].includes(el));
        const hasHydrogen = result.elements.includes('H');
        const hasOxygen = result.elements.includes('O');
        const hasCarbon = result.elements.includes('C');

        // Create smart balancing order that skips already balanced elements initially
        let balancingOrder = [];
        let strategy = "";

        // Check which elements are already balanced with coefficient 1
        const alreadyBalanced = [];
        const needsBalancing = [];

        result.elements.forEach(element => {
            let leftCount = 0;
            let rightCount = 0;

            result.compounds.forEach(compound => {
                const elementCount = compound.elements[element] || 0;
                if (compound.isReactant) {
                    leftCount += elementCount;
                } else {
                    rightCount += elementCount;
                }
            });

            if (leftCount === rightCount) {
                alreadyBalanced.push(element);
            } else {
                needsBalancing.push(element);
            }
        });

        if (hasMetals) {
            strategy = "We balance unbalanced elements first, saving already-balanced elements for when other changes affect them.";
            // Start with unbalanced metals
            balancingOrder = needsBalancing.filter(el => ['Na', 'K', 'Ca', 'Mg', 'Fe', 'Cu', 'Zn', 'Al', 'Ag', 'Pb', 'Sn', 'Cr', 'Mn'].includes(el));
            // Add other unbalanced nonmetals (except H and O)
            balancingOrder = balancingOrder.concat(needsBalancing.filter(el => !['Na', 'K', 'Ca', 'Mg', 'Fe', 'Cu', 'Zn', 'Al', 'Ag', 'Pb', 'Sn', 'Cr', 'Mn', 'H', 'O'].includes(el)));
            // Add hydrogen if unbalanced
            if (needsBalancing.includes('H')) balancingOrder.push('H');
            // Add oxygen if unbalanced  
            if (needsBalancing.includes('O')) balancingOrder.push('O');
            // Finally add any previously balanced elements that might need rebalancing
            balancingOrder = balancingOrder.concat(alreadyBalanced);
        } else if (hasCarbon) {
            strategy = "For organic reactions: balance unbalanced elements first, then recheck others.";
            if (needsBalancing.includes('C')) balancingOrder.push('C');
            if (needsBalancing.includes('H')) balancingOrder.push('H');
            if (needsBalancing.includes('O')) balancingOrder.push('O');
            balancingOrder = balancingOrder.concat(needsBalancing.filter(el => !['C', 'H', 'O'].includes(el)));
            balancingOrder = balancingOrder.concat(alreadyBalanced);
        } else {
            strategy = "We balance unbalanced elements first, then recheck previously balanced ones.";
            balancingOrder = needsBalancing.concat(alreadyBalanced);
        }

        // Simulate the balancing process step by step
        const balancingSteps = this.simulateBalancingProcess(result, balancingOrder);

        return `
            <div class="space-y-6">
                <div class="bg-gradient-to-r from-indigo-50 to-indigo-100 p-4 rounded-lg border-l-4 border-indigo-500">
                    <h4 class="font-bold text-indigo-800 mb-2 flex items-center">
                        <i class="fas fa-route mr-2"></i>Our Balancing Strategy
                    </h4>
                    <p class="text-indigo-700 text-sm mb-3">${strategy}</p>
                    <div class="bg-white p-3 rounded border">
                        <div class="text-sm text-indigo-700">
                            <strong>Watch how the equation changes</strong> as we balance each element in the optimal order!
                        </div>
                    </div>
                </div>

                ${balancingSteps.map((step, index) => `
                    <div class="bg-white border-2 border-indigo-200 rounded-lg overflow-hidden mb-4">
                        <button class="w-full text-left p-4 bg-indigo-100 hover:bg-indigo-150 transition-all duration-300 flex items-center justify-between group" 
                                onclick="const content = this.nextElementSibling; content.style.display = content.style.display === 'none' ? 'block' : 'none'; this.querySelector('.chevron').classList.toggle('rotate-180')">
                            <div class="flex items-center">
                                <div class="w-8 h-8 bg-indigo-600 text-white rounded-full flex items-center justify-center font-bold mr-3 group-hover:scale-110 transition-transform">
                                    ${index + 1}
                                </div>
                                <div>
                                    <h4 class="font-bold text-indigo-800">
                                        ${step.title}
                                    </h4>
                                    <div class="text-sm text-indigo-600 font-medium">
                                        ${step.status}
                                    </div>
                                </div>
                            </div>
                            <i class="fas fa-chevron-down text-indigo-600 chevron transition-transform duration-300"></i>
                        </button>
                        <div class="hidden p-6 bg-white" style="display: none;">
                            <div class="mb-4">
                                <h5 class="font-semibold text-gray-800 mb-2">Before Balancing ${step.element}:</h5>
                                <div class="bg-orange-50 p-4 rounded-lg border-2 border-orange-200">
                                    <div class="font-mono text-2xl text-center font-bold text-gray-800">
                                        ${step.currentEquation}
                                    </div>
                                    <div class="text-center text-sm text-orange-600 mt-2">
                                        ${step.beforeLeftCount !== step.beforeRightCount ?
                `❌ ${step.element}: ${step.beforeLeftCount} ≠ ${step.beforeRightCount} (unbalanced)` :
                `✅ ${step.element}: ${step.beforeLeftCount} = ${step.beforeRightCount} (already balanced)`
            }
                                    </div>
                                </div>
                            </div>

                            <div class="grid grid-cols-1 md:grid-cols-2 gap-6">
                                <div class="bg-blue-50 p-4 rounded-lg border-l-4 border-blue-400">
                                    <h6 class="font-bold text-blue-800 mb-3 flex items-center">
                                        <i class="fas fa-arrow-left mr-2"></i>Left Side Analysis
                                    </h6>
                                    <div class="space-y-2">
                                        ${step.leftCompounds.map(comp => `
                                            <div class="bg-white p-3 rounded border">
                                                <div class="font-mono text-xl font-bold text-blue-700 mb-1">
                                                    ${comp.display}
                                                </div>
                                                <div class="text-sm text-blue-600">
                                                    ${comp.elementBreakdown}
                                                </div>
                                            </div>
                                        `).join('')}
                                    </div>
                                    <div class="mt-3 p-2 bg-blue-100 rounded">
                                        <div class="text-sm font-semibold text-blue-800">
                                            ${step.element}: ${step.leftCount} atoms total
                                        </div>
                                        ${step.leftSideReasoning ? `
                                            <div class="text-xs text-blue-700 mt-1 italic">
                                                ${step.leftSideReasoning}
                                            </div>
                                        ` : ''}
                                    </div>
                                </div>

                                <div class="bg-purple-50 p-4 rounded-lg border-l-4 border-purple-400">
                                    <h6 class="font-bold text-purple-800 mb-3 flex items-center">
                                        <i class="fas fa-arrow-right mr-2"></i>Right Side Analysis
                                    </h6>
                                    <div class="space-y-2">
                                        ${step.rightCompounds.map(comp => `
                                            <div class="bg-white p-3 rounded border">
                                                <div class="font-mono text-xl font-bold text-purple-700 mb-1">
                                                    ${comp.display}
                                                </div>
                                                <div class="text-sm text-purple-600">
                                                    ${comp.elementBreakdown}
                                                </div>
                                            </div>
                                        `).join('')}
                                    </div>
                                    <div class="mt-3 p-2 bg-purple-100 rounded">
                                        <div class="text-sm font-semibold text-purple-800">
                                            ${step.element}: ${step.rightCount} atoms total
                                        </div>
                                        ${step.rightSideReasoning ? `
                                            <div class="text-xs text-purple-700 mt-1 italic">
                                                ${step.rightSideReasoning}
                                            </div>
                                        ` : ''}
                                    </div>
                                </div>
                            </div>

                            <div class="mt-4 p-4 rounded-lg ${step.balanced ? 'bg-green-50 border border-green-200' : 'bg-orange-50 border border-orange-200'}">
                                <div class="flex items-center justify-center">
                                    <div class="text-center">
                                        <div class="text-3xl font-bold ${step.balanced ? 'text-green-600' : 'text-orange-600'} mb-2">
                                            ${step.leftCount} ${step.balanced ? '=' : '≠'} ${step.rightCount}
                                        </div>
                                        ${step.balanced ? `
                                            <div class="mt-4 bg-white p-4 rounded-lg border-2 border-green-300">
                                                <div class="text-sm text-green-700 font-semibold mb-2">After Balancing ${step.element}:</div>
                                                <div class="font-mono text-xl font-bold text-gray-800">
                                                    ${step.progressiveEquation}
                                                </div>
                                                ${step.isLastStep ? `
                                                    <div class="mt-4 bg-gradient-to-r from-green-50 to-green-100 p-4 rounded-lg border-l-4 border-green-500" style="position: relative; z-index: 1;">
                                                        <div class="text-center">
                                                            <div class="w-12 h-12 bg-green-500 text-white rounded-full flex items-center justify-center mx-auto mb-2">
                                                                <i class="fas fa-trophy text-lg"></i>
                                                            </div>
                                                            <h4 class="text-lg font-bold text-green-800 mb-1">🎉 All Elements Balanced!</h4>
                                                            <p class="text-green-700 text-sm">
                                                                Every element now has equal atoms on both sides of the equation!
                                                            </p>
                                                        </div>
                                                    </div>
                                                ` : ''}
                                            </div>
                                        ` : `
                                            <div class="text-sm text-orange-700 font-medium">
                                                ⚖️ Adjusting coefficients to balance ${step.element}...
                                            </div>
                                        `}
                                    </div>
                                </div>
                            </div>

                            ${step.explanation ? `
                                <div class="mt-4 bg-indigo-50 p-4 rounded-lg border border-indigo-200">
                                    <div class="flex items-start">
                                        <i class="fas fa-lightbulb text-indigo-600 mt-1 mr-3"></i>
                                        <div>
                                            <h6 class="font-semibold text-indigo-800 mb-1">What's Happening:</h6>
                                            <p class="text-indigo-700 text-sm">${step.explanation}</p>
                                        </div>
                                    </div>
                                </div>
                            ` : ''}


                        </div>
                    </div>
                `).join('')}
            </div>
        `;
    }

    simulateBalancingProcess(result, balancingOrder) {
        const steps = [];

        // Generate the original unbalanced equation (all coefficients = 1)
        const unbalancedReactants = [];
        const unbalancedProducts = [];
        result.compounds.forEach(compound => {
            if (compound.isReactant) {
                unbalancedReactants.push(compound.formula);
            } else {
                unbalancedProducts.push(compound.formula);
            }
        });
        const unbalancedEquation = `${unbalancedReactants.join(' + ')} → ${unbalancedProducts.join(' + ')}`;

        balancingOrder.forEach((element, index) => {
            // Create "before" coefficients - only apply coefficients for elements balanced BEFORE this step
            const beforeCoeffs = new Array(result.compounds.length).fill(1);
            for (let stepIndex = 0; stepIndex < index; stepIndex++) {
                const balancedElement = balancingOrder[stepIndex];
                result.compounds.forEach((compound, i) => {
                    if (compound.elements[balancedElement]) {
                        beforeCoeffs[i] = result.coefficients[i];
                    }
                });
            }

            // Create "after" coefficients - apply coefficients for elements balanced up to and including this step
            const afterCoeffs = new Array(result.compounds.length).fill(1);
            for (let stepIndex = 0; stepIndex <= index; stepIndex++) {
                const balancedElement = balancingOrder[stepIndex];
                result.compounds.forEach((compound, i) => {
                    if (compound.elements[balancedElement]) {
                        afterCoeffs[i] = result.coefficients[i];
                    }
                });
            }

            // Generate equations for before and after this step
            const beforeReactants = [];
            const beforeProducts = [];
            const afterReactants = [];
            const afterProducts = [];

            result.compounds.forEach((compound, i) => {
                const beforeCoeff = beforeCoeffs[i];
                const afterCoeff = afterCoeffs[i];
                const beforeFormatted = beforeCoeff === 1 ? compound.formula : `${beforeCoeff}${compound.formula}`;
                const afterFormatted = afterCoeff === 1 ? compound.formula : `${afterCoeff}${compound.formula}`;

                if (compound.isReactant) {
                    beforeReactants.push(beforeFormatted);
                    afterReactants.push(afterFormatted);
                } else {
                    beforeProducts.push(beforeFormatted);
                    afterProducts.push(afterFormatted);
                }
            });

            const beforeEquation = `${beforeReactants.join(' + ')} → ${beforeProducts.join(' + ')}`;
            const afterEquation = `${afterReactants.join(' + ')} → ${afterProducts.join(' + ')}`;

            // First calculate atom counts with "before" coefficients to show the imbalance
            let beforeLeftCount = 0;
            let beforeRightCount = 0;

            result.compounds.forEach((compound, i) => {
                const elementCount = compound.elements[element] || 0;
                const coeff = beforeCoeffs[i];
                const totalCount = elementCount * coeff;

                if (compound.isReactant) {
                    beforeLeftCount += totalCount;
                } else {
                    beforeRightCount += totalCount;
                }
            });

            // Calculate atom counts for this element with "after" coefficients
            let leftCount = 0;
            let rightCount = 0;
            let leftCompounds = [];
            let rightCompounds = [];

            result.compounds.forEach((compound, i) => {
                const elementCount = compound.elements[element] || 0;
                const coeff = afterCoeffs[i]; // Use "after" coefficients
                const totalCount = elementCount * coeff;

                const compoundDisplay = coeff === 1 ? compound.formula : `${coeff}${compound.formula}`;
                const elementBreakdown = elementCount > 0 ?
                    `Contains ${elementCount} ${element} atom${elementCount > 1 ? 's' : ''} × ${coeff} = ${totalCount} total` :
                    `No ${element} atoms`;

                if (compound.isReactant) {
                    leftCount += totalCount;
                    if (elementCount > 0) {
                        leftCompounds.push({
                            display: compoundDisplay,
                            elementBreakdown: elementBreakdown
                        });
                    }
                } else {
                    rightCount += totalCount;
                    if (elementCount > 0) {
                        rightCompounds.push({
                            display: compoundDisplay,
                            elementBreakdown: elementBreakdown
                        });
                    }
                }
            });

            // Generate coefficient reasoning for left and right sides
            let leftSideReasoning = "";
            let rightSideReasoning = "";

            if (beforeLeftCount !== beforeRightCount) {
                const imbalance = Math.abs(beforeLeftCount - beforeRightCount);
                const leftChanges = [];
                const rightChanges = [];

                result.compounds.forEach((compound, i) => {
                    if (beforeCoeffs[i] !== afterCoeffs[i]) {
                        const elementCount = compound.elements[element] || 0;
                        if (elementCount > 0) {
                            const coeffChange = afterCoeffs[i] - beforeCoeffs[i];
                            const atomChange = coeffChange * elementCount;
                            const changeInfo = {
                                formula: compound.formula,
                                beforeCoeff: beforeCoeffs[i],
                                afterCoeff: afterCoeffs[i],
                                atomChange: atomChange,
                                elementCount: elementCount
                            };

                            if (compound.isReactant) {
                                leftChanges.push(changeInfo);
                            } else {
                                rightChanges.push(changeInfo);
                            }
                        }
                    }
                });

                if (leftChanges.length > 0) {
                    const totalAtomChange = leftChanges.reduce((sum, change) => sum + change.atomChange, 0);
                    leftSideReasoning = `Need ${totalAtomChange} more ${element} atoms. ` + leftChanges.map(change =>
                        `${change.formula}: ${change.elementCount} ${element} × ${change.afterCoeff} = ${change.elementCount * change.afterCoeff} atoms`
                    ).join(', ');
                } else {
                    leftSideReasoning = `Left side already has correct ${element} count`;
                }

                if (rightChanges.length > 0) {
                    const totalAtomChange = rightChanges.reduce((sum, change) => sum + change.atomChange, 0);
                    rightSideReasoning = `Need ${totalAtomChange} more ${element} atoms. ` + rightChanges.map(change =>
                        `${change.formula}: ${change.elementCount} ${element} × ${change.afterCoeff} = ${change.elementCount * change.afterCoeff} atoms`
                    ).join(', ');
                } else {
                    rightSideReasoning = `Right side already has correct ${element} count`;
                }
            } else {
                // Check if this element was affected by previous balancing steps
                const wasAffected = result.compounds.some((compound, i) =>
                    beforeCoeffs[i] !== 1 && compound.elements[element]
                );

                if (wasAffected) {
                    leftSideReasoning = `${element} rebalanced due to previous coefficient changes`;
                    rightSideReasoning = `${element} rebalanced due to previous coefficient changes`;
                } else {
                    leftSideReasoning = `${element} was already balanced - no changes needed`;
                    rightSideReasoning = `${element} was already balanced - no changes needed`;
                }
            }

            // Determine explanation based on element type and position
            let explanation = "";
            let title = "";
            let status = "";

            if (['Na', 'K', 'Ca', 'Mg', 'Fe', 'Cu', 'Zn', 'Al', 'Ag', 'Pb', 'Sn', 'Cr', 'Mn'].includes(element)) {
                title = `Balancing ${element} (Metal)`;
                status = "Step " + (index + 1) + " of " + balancingOrder.length;
                explanation = `Metals like ${element} usually appear in fewer compounds, making them easier to balance first. We adjust coefficients to make ${element} atoms equal on both sides.`;
            } else if (element === 'C') {
                title = `Balancing ${element} (Carbon Backbone)`;
                status = "Step " + (index + 1) + " of " + balancingOrder.length;
                explanation = `Carbon forms the backbone of organic molecules. We balance it first to establish the correct molecular ratios for the entire reaction.`;
            } else if (element === 'H') {
                title = `Balancing ${element} (Hydrogen)`;
                status = "Step " + (index + 1) + " of " + balancingOrder.length;
                explanation = `Hydrogen appears in many compounds (water, acids, bases). We balance it after other elements are set, but before oxygen.`;
            } else if (element === 'O') {
                title = `Balancing ${element} (Oxygen - Final Step)`;
                status = "Final Step";
                explanation = `Oxygen is saved for last because it appears in the most compounds. By now, the coefficients we've set for other elements usually balance oxygen automatically!`;
            } else {
                title = `Balancing ${element} (Nonmetal)`;
                status = "Step " + (index + 1) + " of " + balancingOrder.length;
                explanation = `This nonmetal is balanced after metals but before hydrogen and oxygen, following our systematic approach.`;
            }

            steps.push({
                element: element,
                title: title,
                status: status,
                currentEquation: beforeEquation, // Show equation at start of this step
                leftCompounds: leftCompounds, // These show "after" coefficients in the analysis
                rightCompounds: rightCompounds, // These show "after" coefficients in the analysis
                leftCount: leftCount, // "After" atom counts
                rightCount: rightCount, // "After" atom counts
                beforeLeftCount: beforeLeftCount, // "Before" atom counts for display
                beforeRightCount: beforeRightCount, // "Before" atom counts for display
                balanced: leftCount === rightCount, // Should be true when this element is balanced
                explanation: explanation,
                leftSideReasoning: leftSideReasoning, // Left side coefficient reasoning
                rightSideReasoning: rightSideReasoning, // Right side coefficient reasoning
                isLastStep: index === balancingOrder.length - 1,
                progressiveEquation: afterEquation // Show equation after this step is balanced
            });
        });

        return steps;
    }

    generateApplyStep(result) {
        return `
            <div class="space-y-4">
                <div class="bg-gradient-to-r from-orange-50 to-orange-100 p-4 rounded-lg border-l-4 border-orange-500">
                    <h4 class="font-bold text-orange-800 mb-2 flex items-center">
                        <i class="fas fa-calculator mr-2"></i>Applying Our Coefficients
                    </h4>
                    <p class="text-orange-700 text-sm">
                        Now we take the coefficients we calculated and put them in front of each compound.
                    </p>
                </div>
                
                <div class="bg-white border-2 border-orange-200 rounded-lg p-4">
                    <h4 class="font-bold text-orange-800 mb-4 flex items-center">
                        <i class="fas fa-arrow-right mr-2"></i>Before → After
                    </h4>
                    
                    <div class="space-y-4">
                        <div class="bg-slate-50 p-4 rounded-lg border-l-4 border-slate-400">
                            <h5 class="font-bold text-slate-700 mb-2">BEFORE (Unbalanced):</h5>
                            <div class="font-mono text-lg text-center bg-white p-3 rounded border">
                                ${result.compounds.filter(c => c.isReactant).map(c => c.formula).join(' + ')} → ${result.compounds.filter(c => !c.isReactant).map(c => c.formula).join(' + ')}
                            </div>
                        </div>
                        
                        <div class="flex justify-center">
                            <div class="w-12 h-12 bg-orange-500 text-white rounded-full flex items-center justify-center">
                                <i class="fas fa-arrow-down text-xl"></i>
                            </div>
                        </div>
                        
                        <div class="bg-emerald-50 p-4 rounded-lg border-l-4 border-emerald-400">
                            <h5 class="font-bold text-emerald-700 mb-2">AFTER (Balanced):</h5>
                            <div class="font-mono text-xl text-center bg-white p-4 rounded border font-bold">
                                ${result.balancedEquation}
                            </div>
                        </div>
                    </div>
                </div>
                
                <div class="bg-white border-2 border-orange-200 rounded-lg p-4">
                    <h4 class="font-bold text-orange-800 mb-4 flex items-center">
                        <i class="fas fa-trophy mr-2"></i>The Coefficients We Found
                    </h4>
                    <div class="grid gap-3">
                        ${result.compounds.map((compound, i) => {
            const isReactant = compound.isReactant;
            const bgColor = isReactant ? 'bg-red-50 border-red-200' : 'bg-green-50 border-green-200';
            const textColor = isReactant ? 'text-red-700' : 'text-green-700';
            const coefficient = result.coefficients[i];

            return `
                                <div class="${bgColor} border-2 rounded-lg p-4 flex items-center justify-between">
                                    <div class="flex items-center">
                                        <div class="w-12 h-12 bg-orange-500 text-white rounded-full flex items-center justify-center font-bold text-xl mr-4">
                                            ${coefficient}
                                        </div>
                                        <span class="font-mono text-xl font-bold ${textColor}">${compound.formula}</span>
                                    </div>
                                    <div class="text-sm ${textColor}">
                                        ${coefficient === 1 ? 'No coefficient needed (1 is implied)' : `We need ${coefficient} molecules`}
                                    </div>
                                </div>
                            `;
        }).join('')}
                    </div>
                </div>
                
                <div class="bg-orange-50 p-4 rounded-lg border border-orange-200">
                    <div class="flex items-start">
                        <i class="fas fa-lightbulb text-orange-600 mt-1 mr-3"></i>
                        <div>
                            <h5 class="font-semibold text-orange-800 mb-1">What These Numbers Mean</h5>
                            <p class="text-orange-700 text-sm">
                                Each coefficient tells us how many molecules (or moles) of that compound we need. 
                                For example, "2H₂O" means we need 2 water molecules for this reaction to work perfectly.
                            </p>
                        </div>
                    </div>
                </div>
            </div>
        `;
    }

    generateMatrixStep(result) {
        return `
            <div class="space-y-4">
                <div class="bg-gradient-to-r from-indigo-50 to-indigo-100 p-4 rounded-lg border-l-4 border-indigo-500">
                    <h4 class="font-bold text-indigo-800 mb-2 flex items-center">
                        <i class="fas fa-balance-scale mr-2"></i>The Balance Equations
                    </h4>
                    <p class="text-indigo-700 text-sm">
                        For each element, we create an equation where: (Reactant atoms) - (Product atoms) = 0
                    </p>
                </div>
                
                <div class="bg-white border-2 border-indigo-200 rounded-lg overflow-hidden">
                    <div class="bg-indigo-100 p-3 border-b border-indigo-200">
                        <h4 class="font-bold text-indigo-800 flex items-center">
                            <i class="fas fa-table mr-2"></i>Balance Matrix
                        </h4>
                    </div>
                    <div class="overflow-x-auto">
                        <table class="w-full">
                            <thead class="bg-gray-50">
                                <tr>
                                    <th class="text-left p-3 font-bold text-gray-700 border-r">Element</th>
                                    ${result.compounds.map((compound, i) => {
            const isReactant = compound.isReactant;
            const bgColor = isReactant ? 'bg-red-100' : 'bg-green-100';
            const textColor = isReactant ? 'text-red-700' : 'text-green-700';
            return `<th class="text-center p-3 font-bold ${textColor} ${bgColor} border-r">${compound.formula}</th>`;
        }).join('')}
                                    <th class="text-center p-3 font-bold text-gray-700">=</th>
                                </tr>
                            </thead>
                            <tbody>
                                ${result.elements.map((element, i) => {
            return `
                                        <tr class="border-t hover:bg-gray-50">
                                            <td class="p-3 font-bold text-gray-800 border-r bg-gray-50">${element}</td>
                                            ${result.matrix[i].map((coeff, j) => {
                const isReactant = result.compounds[j].isReactant;
                const sign = coeff > 0 ? '+' : '';
                const color = coeff > 0 ? 'text-red-600' : 'text-green-600';
                const bgColor = isReactant ? 'bg-red-50' : 'bg-green-50';
                return `<td class="text-center p-3 font-mono font-bold ${color} ${bgColor} border-r">${sign}${coeff}</td>`;
            }).join('')}
                                            <td class="text-center p-3 font-bold text-gray-800">0</td>
                                        </tr>
                                    `;
        }).join('')}
                            </tbody>
                        </table>
                    </div>
                </div>
                
                <div class="bg-indigo-50 p-4 rounded-lg border border-indigo-200">
                    <div class="flex items-start">
                        <i class="fas fa-lightbulb text-indigo-600 mt-1 mr-3"></i>
                        <div>
                            <h5 class="font-semibold text-indigo-800 mb-1">Reading the Matrix</h5>
                            <p class="text-indigo-700 text-sm mb-2">
                                Each row represents one element. Positive numbers are reactants, negative numbers are products.
                                Each equation must equal zero for the reaction to be balanced.
                            </p>
                            <div class="text-xs text-indigo-600">
                                <span class="inline-block bg-red-100 text-red-700 px-2 py-1 rounded mr-2">+ Reactants</span>
                                <span class="inline-block bg-green-100 text-green-700 px-2 py-1 rounded">- Products</span>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        `;
    }

    generateSolutionStep(result) {
        return `
            <div class="space-y-4">
                <div class="bg-gradient-to-r from-orange-50 to-orange-100 p-4 rounded-lg border-l-4 border-orange-500">
                    <h4 class="font-bold text-orange-800 mb-2 flex items-center">
                        <i class="fas fa-calculator mr-2"></i>Mathematical Solution
                    </h4>
                    <p class="text-orange-700 text-sm">
                        Using linear algebra (Gaussian elimination), we solve the system of equations to find the smallest whole number coefficients.
                    </p>
                </div>
                
                <div class="bg-white border-2 border-orange-200 rounded-lg p-4">
                    <h4 class="font-bold text-orange-800 mb-4 flex items-center">
                        <i class="fas fa-trophy mr-2"></i>Calculated Coefficients
                    </h4>
                    <div class="grid gap-3">
                        ${result.compounds.map((compound, i) => {
            const isReactant = compound.isReactant;
            const bgColor = isReactant ? 'bg-red-50 border-red-200' : 'bg-green-50 border-green-200';
            const textColor = isReactant ? 'text-red-700' : 'text-green-700';
            const icon = isReactant ? 'fas fa-arrow-right' : 'fas fa-bullseye';
            const coefficient = result.coefficients[i];

            return `
                                <div class="${bgColor} border-2 rounded-lg p-4 flex items-center justify-between">
                                    <div class="flex items-center">
                                        <i class="${icon} ${textColor} mr-3"></i>
                                        <span class="font-mono text-xl font-bold ${textColor}">${compound.formula}</span>
                                    </div>
                                    <div class="flex items-center">
                                        <span class="text-sm text-gray-600 mr-2">coefficient =</span>
                                        <div class="w-12 h-12 bg-orange-500 text-white rounded-full flex items-center justify-center font-bold text-xl">
                                            ${coefficient}
                                        </div>
                                    </div>
                                </div>
                            `;
        }).join('')}
                    </div>
                </div>
                
                <div class="bg-orange-100 border-2 border-orange-300 rounded-lg p-4">
                    <h4 class="font-bold text-orange-800 mb-3 flex items-center">
                        <i class="fas fa-magic mr-2"></i>Applying the Coefficients
                    </h4>
                    <div class="bg-white p-4 rounded-lg border-2 border-orange-200">
                        <div class="text-center">
                            <div class="text-2xl font-mono font-bold text-gray-800 mb-2">
                                ${result.compounds.filter(c => c.isReactant).map((compound, i) => {
            const originalIndex = result.compounds.indexOf(compound);
            const coeff = result.coefficients[originalIndex];
            const formatted = coeff === 1 ? compound.formula : `${coeff}${compound.formula}`;
            return `<span class="text-red-600">${formatted}</span>`;
        }).join(' + ')} → ${result.compounds.filter(c => !c.isReactant).map((compound, i) => {
            const originalIndex = result.compounds.indexOf(compound);
            const coeff = result.coefficients[originalIndex];
            const formatted = coeff === 1 ? compound.formula : `${coeff}${compound.formula}`;
            return `<span class="text-green-600">${formatted}</span>`;
        }).join(' + ')}
                            </div>
                            <p class="text-sm text-gray-600">Balanced Chemical Equation</p>
                        </div>
                    </div>
                </div>
                
                <div class="bg-orange-50 p-4 rounded-lg border border-orange-200">
                    <div class="flex items-start">
                        <i class="fas fa-lightbulb text-orange-600 mt-1 mr-3"></i>
                        <div>
                            <h5 class="font-semibold text-orange-800 mb-1">How we got these numbers</h5>
                            <p class="text-orange-700 text-sm">
                                The algorithm uses Gaussian elimination to solve the matrix equations, then converts the results 
                                to the smallest possible whole numbers. This ensures we use the minimum amount of each compound.
                            </p>
                        </div>
                    </div>
                </div>
            </div>
        `;
    }

    generateVerificationStep(result) {
        const elementColors = {
            'H': 'bg-red-100 text-red-800 border-red-300',
            'O': 'bg-blue-100 text-blue-800 border-blue-300',
            'C': 'bg-gray-100 text-gray-800 border-gray-300',
            'N': 'bg-green-100 text-green-800 border-green-300',
            'S': 'bg-yellow-100 text-yellow-800 border-yellow-300',
            'Cl': 'bg-purple-100 text-purple-800 border-purple-300',
            'Na': 'bg-orange-100 text-orange-800 border-orange-300',
            'Ca': 'bg-pink-100 text-pink-800 border-pink-300',
            'Fe': 'bg-indigo-100 text-indigo-800 border-indigo-300'
        };

        return `
            <div class="space-y-4">
                <div class="bg-gradient-to-r from-green-50 to-green-100 p-4 rounded-lg border-l-4 border-green-500">
                    <h4 class="font-bold text-green-800 mb-2 flex items-center">
                        <i class="fas fa-balance-scale mr-2"></i>Final Verification
                    </h4>
                    <p class="text-green-700 text-sm">
                        Let's double-check that each element has the same number of atoms on both sides of the equation.
                    </p>
                </div>
                
                <div class="bg-white border-2 border-green-200 rounded-lg p-4">
                    <h4 class="font-bold text-green-800 mb-4 flex items-center">
                        <i class="fas fa-microscope mr-2"></i>Atom Count Verification
                    </h4>
                    <div class="space-y-3">
                        ${result.elements.map(element => {
            let reactantCount = 0;
            let productCount = 0;

            result.compounds.forEach((compound, i) => {
                const count = (compound.elements[element] || 0) * result.coefficients[i];
                if (compound.isReactant) {
                    reactantCount += count;
                } else {
                    productCount += count;
                }
            });

            const colorClass = elementColors[element] || 'bg-gray-100 text-gray-800 border-gray-300';
            const isBalanced = reactantCount === productCount;

            return `
                                <div class="bg-gray-50 border-2 border-gray-200 rounded-lg p-4">
                                    <div class="flex items-center justify-between">
                                        <div class="flex items-center">
                                            <span class="w-8 h-8 rounded-full border-2 ${colorClass} flex items-center justify-center text-sm font-bold mr-3">
                                                ${element}
                                            </span>
                                            <span class="font-bold text-gray-800">${element} atoms:</span>
                                        </div>
                                        <div class="flex items-center space-x-4">
                                            <div class="text-center">
                                                <div class="bg-red-100 text-red-800 px-3 py-1 rounded-full font-bold">
                                                    ${reactantCount}
                                                </div>
                                                <div class="text-xs text-red-600 mt-1">Reactants</div>
                                            </div>
                                            <div class="text-2xl font-bold text-gray-400">=</div>
                                            <div class="text-center">
                                                <div class="bg-green-100 text-green-800 px-3 py-1 rounded-full font-bold">
                                                    ${productCount}
                                                </div>
                                                <div class="text-xs text-green-600 mt-1">Products</div>
                                            </div>
                                            <div class="text-2xl ${isBalanced ? 'text-green-500' : 'text-red-500'}">
                                                <i class="fas fa-${isBalanced ? 'check-circle' : 'times-circle'}"></i>
                                            </div>
                                        </div>
                                    </div>
                                </div>
                            `;
        }).join('')}
                    </div>
                </div>
                
                <div class="bg-gradient-to-r from-green-100 to-green-200 border-2 border-green-300 rounded-lg p-6">
                    <div class="text-center">
                        <div class="w-16 h-16 bg-green-500 text-white rounded-full flex items-center justify-center mx-auto mb-4">
                            <i class="fas fa-trophy text-2xl"></i>
                        </div>
                        <h3 class="text-2xl font-bold text-green-800 mb-2">🎉 Equation Successfully Balanced!</h3>
                        <div class="bg-white p-4 rounded-lg border-2 border-green-300 mb-4">
                            <div class="text-3xl font-mono font-bold text-gray-800">
                                ${result.balancedEquation}
                            </div>
                        </div>
                        <p class="text-green-700 font-medium">
                            All elements are now balanced according to the Law of Conservation of Mass!
                        </p>
                    </div>
                </div>
                
                <div class="bg-green-50 p-4 rounded-lg border border-green-200">
                    <div class="flex items-start">
                        <i class="fas fa-lightbulb text-green-600 mt-1 mr-3"></i>
                        <div>
                            <h5 class="font-semibold text-green-800 mb-1">Why this matters</h5>
                            <p class="text-green-700 text-sm">
                                A balanced equation follows the Law of Conservation of Mass - matter cannot be created or destroyed 
                                in a chemical reaction. The same number of each type of atom must appear on both sides.
                            </p>
                        </div>
                    </div>
                </div>
            </div>
        `;
    }

    showError(message) {
        console.log('Showing error:', message);
        const errorSection = document.getElementById('error-section');
        const errorMessage = document.getElementById('error-message');

        if (errorSection && errorMessage) {
            errorMessage.textContent = message;
            errorSection.classList.remove('hidden');
            errorSection.scrollIntoView({ behavior: 'smooth' });
        }
    }

    hideError() {
        const errorSection = document.getElementById('error-section');
        if (errorSection) {
            errorSection.classList.add('hidden');
        }
    }

    hideResults() {
        const resultsSection = document.getElementById('results-section');
        if (resultsSection) {
            resultsSection.classList.add('hidden');
        }
    }

    clearInputs() {
        const reactantsInput = document.getElementById('reactants-input');
        const productsInput = document.getElementById('products-input');

        if (reactantsInput) {
            reactantsInput.value = '';
            reactantsInput.focus();
        }
        if (productsInput) {
            productsInput.value = '';
        }

        // Reset dropdown filters
        this.filterDropdown('reactants', '');
        this.filterDropdown('products', '');

        this.closeAllDropdowns();
        this.hideError();
        this.hideResults();
    }
}

// Clean Dark Mode Implementation
const toggleButton = document.getElementById('theme-toggle');
const themeIcon = document.getElementById('theme-icon');
const body = document.body;

// Load saved theme
if (localStorage.getItem('theme') === 'dark') {
    body.classList.add('dark-mode');
    if (themeIcon) themeIcon.className = 'fas fa-sun';
}

// When button is clicked, toggle theme
if (toggleButton) {
    toggleButton.addEventListener('click', () => {
        body.classList.toggle('dark-mode');

        // Update icon
        const isDark = body.classList.contains('dark-mode');
        if (themeIcon) {
            themeIcon.className = isDark ? 'fas fa-sun' : 'fas fa-moon';
        }

        // Save theme choice
        if (body.classList.contains('dark-mode')) {
            localStorage.setItem('theme', 'dark');
        } else {
            localStorage.setItem('theme', 'light');
        }
    });
}

// Initialize the application
document.addEventListener('DOMContentLoaded', () => {
    console.log('DOM loaded, initializing app...');
    new EquationBalancerUI();
});