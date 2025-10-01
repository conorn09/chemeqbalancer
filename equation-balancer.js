class ChemicalEquationBalancer {
    constructor() {
        this.elements = [];
        this.compounds = [];
        this.matrix = [];
        this.solution = [];
    }

    // Parse chemical formula to extract elements and their counts
    parseFormula(formula) {
        const elements = {};

        // Convert subscript numbers to regular numbers
        const normalizedFormula = formula
            .replace(/₀/g, '0').replace(/₁/g, '1').replace(/₂/g, '2')
            .replace(/₃/g, '3').replace(/₄/g, '4').replace(/₅/g, '5')
            .replace(/₆/g, '6').replace(/₇/g, '7').replace(/₈/g, '8')
            .replace(/₉/g, '9');

        const regex = /([A-Z][a-z]?)(\d*)/g;
        let match;

        while ((match = regex.exec(normalizedFormula)) !== null) {
            const element = match[1];
            const count = match[2] ? parseInt(match[2]) : 1;
            elements[element] = (elements[element] || 0) + count;
        }

        return elements;
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

    // Solve the system using Gaussian elimination
    solveMatrix() {
        const matrix = this.matrix.map(row => [...row]);
        const numRows = matrix.length;
        const numCols = matrix[0].length;

        // Gaussian elimination
        for (let i = 0; i < Math.min(numRows, numCols - 1); i++) {
            // Find pivot
            let maxRow = i;
            for (let k = i + 1; k < numRows; k++) {
                if (Math.abs(matrix[k][i]) > Math.abs(matrix[maxRow][i])) {
                    maxRow = k;
                }
            }

            // Swap rows
            [matrix[i], matrix[maxRow]] = [matrix[maxRow], matrix[i]];

            // Make all rows below this one 0 in current column
            for (let k = i + 1; k < numRows; k++) {
                if (matrix[i][i] !== 0) {
                    const factor = matrix[k][i] / matrix[i][i];
                    for (let j = i; j < numCols; j++) {
                        matrix[k][j] -= factor * matrix[i][j];
                    }
                }
            }
        }

        // Back substitution to find a particular solution
        const solution = new Array(numCols).fill(0);
        solution[numCols - 1] = 1; // Set last coefficient to 1

        for (let i = Math.min(numRows, numCols - 1) - 1; i >= 0; i--) {
            let sum = 0;
            for (let j = i + 1; j < numCols; j++) {
                sum += matrix[i][j] * solution[j];
            }
            if (matrix[i][i] !== 0) {
                solution[i] = -sum / matrix[i][i];
            }
        }

        // Convert to positive integers
        const gcd = (a, b) => b === 0 ? a : gcd(b, a % b);
        const lcm = (a, b) => Math.abs(a * b) / gcd(a, b);

        // Find LCM of denominators to make all coefficients integers
        let denominator = 1;
        const fractions = solution.map(x => {
            if (x === 0) return { num: 0, den: 1 };
            const precision = 1000000;
            const num = Math.round(x * precision);
            const den = precision;
            const g = gcd(Math.abs(num), den);
            return { num: num / g, den: den / g };
        });

        fractions.forEach(frac => {
            if (frac.den !== 1) {
                denominator = lcm(denominator, frac.den);
            }
        });

        const integerSolution = fractions.map(frac =>
            Math.round((frac.num * denominator) / frac.den)
        );

        // Make all coefficients positive
        const minCoeff = Math.min(...integerSolution.filter(x => x !== 0));
        if (minCoeff < 0) {
            integerSolution.forEach((coeff, i) => {
                integerSolution[i] = -coeff;
            });
        }

        // Reduce to smallest integers
        const coeffGcd = integerSolution.reduce((acc, val) =>
            val === 0 ? acc : gcd(acc, Math.abs(val)), 0
        );

        if (coeffGcd > 1) {
            integerSolution.forEach((coeff, i) => {
                integerSolution[i] = Math.round(coeff / coeffGcd);
            });
        }

        this.solution = integerSolution;
        return integerSolution;
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
            let sum = 0;
            for (let j = 0; j < this.compounds.length; j++) {
                const element = this.elements[i];
                const count = this.compounds[j].elements[element] || 0;
                const coeff = coefficients[j];
                sum += (this.compounds[j].isReactant ? count : -count) * coeff;
            }
            if (Math.abs(sum) > 1e-10) {
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
    { formula: "SO₂", name: "Sulfur dioxide", category: "Gas" },
    { formula: "H₂S", name: "Hydrogen sulfide", category: "Gas" },
    { formula: "NO", name: "Nitric oxide", category: "Gas" },
    { formula: "NO₂", name: "Nitrogen dioxide", category: "Gas" },
    { formula: "Cl₂", name: "Chlorine gas", category: "Gas" },

    // Water and common liquids
    { formula: "H₂O", name: "Water", category: "Liquid" },
    { formula: "H₂O₂", name: "Hydrogen peroxide", category: "Liquid" },
    { formula: "C₂H₅OH", name: "Ethanol", category: "Liquid" },
    { formula: "CH₃OH", name: "Methanol", category: "Liquid" },

    // Common acids
    { formula: "HCl", name: "Hydrochloric acid", category: "Acid" },
    { formula: "H₂SO₄", name: "Sulfuric acid", category: "Acid" },
    { formula: "HNO₃", name: "Nitric acid", category: "Acid" },
    { formula: "CH₃COOH", name: "Acetic acid", category: "Acid" },
    { formula: "H₃PO₄", name: "Phosphoric acid", category: "Acid" },
    { formula: "HF", name: "Hydrofluoric acid", category: "Acid" },
    { formula: "HBr", name: "Hydrobromic acid", category: "Acid" },
    { formula: "HI", name: "Hydroiodic acid", category: "Acid" },

    // Common bases
    { formula: "NaOH", name: "Sodium hydroxide", category: "Base" },
    { formula: "KOH", name: "Potassium hydroxide", category: "Base" },
    { formula: "Ca(OH)₂", name: "Calcium hydroxide", category: "Base" },
    { formula: "Mg(OH)₂", name: "Magnesium hydroxide", category: "Base" },
    { formula: "Ba(OH)₂", name: "Barium hydroxide", category: "Base" },

    // Common salts
    { formula: "NaCl", name: "Sodium chloride", category: "Salt" },
    { formula: "KCl", name: "Potassium chloride", category: "Salt" },
    { formula: "CaCl₂", name: "Calcium chloride", category: "Salt" },
    { formula: "MgCl₂", name: "Magnesium chloride", category: "Salt" },
    { formula: "Na₂SO₄", name: "Sodium sulfate", category: "Salt" },
    { formula: "K₂SO₄", name: "Potassium sulfate", category: "Salt" },
    { formula: "CaSO₄", name: "Calcium sulfate", category: "Salt" },
    { formula: "Na₂CO₃", name: "Sodium carbonate", category: "Salt" },
    { formula: "CaCO₃", name: "Calcium carbonate", category: "Salt" },
    { formula: "NaHCO₃", name: "Sodium bicarbonate", category: "Salt" },
    { formula: "AgNO₃", name: "Silver nitrate", category: "Salt" },
    { formula: "Pb(NO₃)₂", name: "Lead nitrate", category: "Salt" },

    // Common oxides
    { formula: "Fe₂O₃", name: "Iron(III) oxide", category: "Oxide" },
    { formula: "FeO", name: "Iron(II) oxide", category: "Oxide" },
    { formula: "Al₂O₃", name: "Aluminum oxide", category: "Oxide" },
    { formula: "CuO", name: "Copper(II) oxide", category: "Oxide" },
    { formula: "Cu₂O", name: "Copper(I) oxide", category: "Oxide" },
    { formula: "ZnO", name: "Zinc oxide", category: "Oxide" },
    { formula: "MgO", name: "Magnesium oxide", category: "Oxide" },
    { formula: "CaO", name: "Calcium oxide", category: "Oxide" },
    { formula: "Na₂O", name: "Sodium oxide", category: "Oxide" },
    { formula: "K₂O", name: "Potassium oxide", category: "Oxide" },

    // Elements
    { formula: "Fe", name: "Iron", category: "Element" },
    { formula: "Cu", name: "Copper", category: "Element" },
    { formula: "Zn", name: "Zinc", category: "Element" },
    { formula: "Al", name: "Aluminum", category: "Element" },
    { formula: "Mg", name: "Magnesium", category: "Element" },
    { formula: "Ca", name: "Calcium", category: "Element" },
    { formula: "Na", name: "Sodium", category: "Element" },
    { formula: "K", name: "Potassium", category: "Element" },
    { formula: "Ag", name: "Silver", category: "Element" },
    { formula: "Pb", name: "Lead", category: "Element" },
    { formula: "C", name: "Carbon", category: "Element" },
    { formula: "S", name: "Sulfur", category: "Element" },
    { formula: "P", name: "Phosphorus", category: "Element" }
];

// UI Controller
class EquationBalancerUI {
    constructor() {
        this.balancer = new ChemicalEquationBalancer();
        this.compounds = COMPOUND_DATABASE;
        this.activeDropdown = null;
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
                this.showDropdown('reactants');
                // Filter based on current input
                const query = this.getCurrentCompound(reactantsInput.value);
                this.filterDropdown('reactants', query);
            });

            reactantsInput.addEventListener('input', (e) => {
                const query = this.getCurrentCompound(e.target.value);
                this.showDropdown('reactants');
                this.filterDropdown('reactants', query);
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
                this.showDropdown('products');
                // Filter based on current input
                const query = this.getCurrentCompound(productsInput.value);
                this.filterDropdown('products', query);
            });

            productsInput.addEventListener('input', (e) => {
                const query = this.getCurrentCompound(e.target.value);
                this.showDropdown('products');
                this.filterDropdown('products', query);
            });
        }
    }

    // Get the current compound being typed (after the last +)
    getCurrentCompound(inputValue) {
        return inputValue.split('+').pop().trim();
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
        const compounds = currentValue.split('+').map(c => c.trim()).filter(c => c.length > 0);

        // Replace the last compound (the one being typed) with the selected formula
        if (compounds.length === 0 || currentValue === '') {
            input.value = formula;
        } else {
            // Check if we're replacing an incomplete compound or adding a new one
            const lastCompound = compounds[compounds.length - 1];
            const currentTyping = this.getCurrentCompound(currentValue);

            if (currentTyping && currentTyping.length > 0) {
                // Replace the compound being typed
                compounds[compounds.length - 1] = formula;
            } else {
                // Add new compound
                compounds.push(formula);
            }

            input.value = compounds.join(' + ');
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
        const dropdown = document.getElementById(`${type}-dropdown`);
        if (dropdown) {
            dropdown.classList.remove('hidden');
            this.activeDropdown = type;
        }
    }

    closeDropdown(type) {
        const dropdown = document.getElementById(`${type}-dropdown`);
        if (dropdown) {
            dropdown.classList.add('hidden');
            if (this.activeDropdown === type) {
                this.activeDropdown = null;
            }
        }
    }

    closeAllDropdowns() {
        this.closeDropdown('reactants');
        this.closeDropdown('products');
    }

    addExamplePlaceholders() {
        const examples = [
            { reactants: 'H2 + O2', products: 'H2O' },
            { reactants: 'CH4 + O2', products: 'CO2 + H2O' },
            { reactants: 'C2H6 + O2', products: 'CO2 + H2O' },
            { reactants: 'Fe + O2', products: 'Fe2O3' },
            { reactants: 'NH3 + O2', products: 'NO + H2O' }
        ];

        const reactantsInput = document.getElementById('reactants-input');
        const productsInput = document.getElementById('products-input');

        if (reactantsInput && productsInput) {
            const currentExample = examples[Math.floor(Math.random() * examples.length)];

            reactantsInput.addEventListener('focus', function () {
                if (!this.value) {
                    this.placeholder = currentExample.reactants;
                }
            });

            productsInput.addEventListener('focus', function () {
                if (!this.value) {
                    this.placeholder = currentExample.products;
                }
            });
        }
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
                title: '🔍 Step 1: Identify What We\'re Balancing',
                icon: 'fas fa-search',
                color: 'blue',
                content: this.generateIdentifyStep(originalEquation, result)
            },
            {
                title: '🎯 Step 2: Count Atoms (Before Balancing)',
                icon: 'fas fa-atom',
                color: 'purple',
                content: this.generateCountStep(result)
            },
            {
                title: '⚖️ Step 3: Balance Step by Step',
                icon: 'fas fa-balance-scale',
                color: 'indigo',
                content: this.generateBalancingStep(result, originalEquation)
            },
            {
                title: '🧮 Step 4: Apply the Coefficients',
                icon: 'fas fa-calculator',
                color: 'orange',
                content: this.generateApplyStep(result)
            },
            {
                title: '✅ Step 5: Double-Check Our Work',
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
                <div class="p-6 bg-white border-t-2 border-${step.color}-100 ${index === 0 ? '' : 'hidden'}">
                    ${step.content}
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
        const hasMetals = result.elements.some(el => ['Na', 'K', 'Ca', 'Mg', 'Fe', 'Cu', 'Zn', 'Al', 'Ag', 'Pb'].includes(el));
        const hasHydrogen = result.elements.includes('H');
        const hasOxygen = result.elements.includes('O');
        const hasCarbon = result.elements.includes('C');

        // Create balancing order based on best practices
        let balancingOrder = [];
        let strategy = "";

        if (hasMetals) {
            strategy = "We'll start with metals (they're usually easiest), then nonmetals, then hydrogen, and oxygen last.";
            balancingOrder = result.elements.filter(el => ['Na', 'K', 'Ca', 'Mg', 'Fe', 'Cu', 'Zn', 'Al', 'Ag', 'Pb'].includes(el));
            balancingOrder = balancingOrder.concat(result.elements.filter(el => !['Na', 'K', 'Ca', 'Mg', 'Fe', 'Cu', 'Zn', 'Al', 'Ag', 'Pb', 'H', 'O'].includes(el)));
            if (hasHydrogen) balancingOrder.push('H');
            if (hasOxygen) balancingOrder.push('O');
        } else if (hasCarbon) {
            strategy = "This looks like an organic reaction. We'll balance carbon first, then hydrogen, then oxygen.";
            balancingOrder = ['C'];
            if (hasHydrogen) balancingOrder.push('H');
            if (hasOxygen) balancingOrder.push('O');
            balancingOrder = balancingOrder.concat(result.elements.filter(el => !['C', 'H', 'O'].includes(el)));
        } else {
            strategy = "We'll balance elements that appear in fewer compounds first, saving diatomic molecules (like O₂) for last.";
            balancingOrder = [...result.elements];
        }

        return `
            <div class="space-y-4">
                <div class="bg-gradient-to-r from-indigo-50 to-indigo-100 p-4 rounded-lg border-l-4 border-indigo-500">
                    <h4 class="font-bold text-indigo-800 mb-2 flex items-center">
                        <i class="fas fa-strategy mr-2"></i>Our Balancing Strategy
                    </h4>
                    <p class="text-indigo-700 text-sm">${strategy}</p>
                </div>
                
                <div class="bg-white border-2 border-indigo-200 rounded-lg p-4">
                    <h4 class="font-bold text-indigo-800 mb-4 flex items-center">
                        <i class="fas fa-list-ol mr-2"></i>Step-by-Step Balancing Order
                    </h4>
                    <div class="space-y-3">
                        ${balancingOrder.map((element, index) => {
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
            const colorClass = elementColors[element] || 'bg-gray-100 text-gray-800 border-gray-300';

            let explanation = "";
            if (['Na', 'K', 'Ca', 'Mg', 'Fe', 'Cu', 'Zn', 'Al', 'Ag', 'Pb'].includes(element)) {
                explanation = "Metal - usually appears in fewer compounds, easier to balance first";
            } else if (element === 'C') {
                explanation = "Carbon - backbone of organic molecules, balance first";
            } else if (element === 'H') {
                explanation = "Hydrogen - balance after other elements but before oxygen";
            } else if (element === 'O') {
                explanation = "Oxygen - save for last, often appears in many compounds";
            } else {
                explanation = "Nonmetal - balance after metals but before H and O";
            }

            return `
                                <div class="flex items-center p-3 bg-gray-50 rounded-lg border">
                                    <div class="w-8 h-8 bg-indigo-500 text-white rounded-full flex items-center justify-center font-bold mr-4">
                                        ${index + 1}
                                    </div>
                                    <span class="w-10 h-10 rounded-full border-2 ${colorClass} flex items-center justify-center text-sm font-bold mr-4">
                                        ${element}
                                    </span>
                                    <div>
                                        <div class="font-semibold text-gray-800">${element}</div>
                                        <div class="text-sm text-gray-600">${explanation}</div>
                                    </div>
                                </div>
                            `;
        }).join('')}
                    </div>
                </div>
                
                <div class="bg-indigo-100 border-2 border-indigo-300 rounded-lg p-4">
                    <h4 class="font-bold text-indigo-800 mb-3 flex items-center">
                        <i class="fas fa-magic mr-2"></i>The Math Behind It
                    </h4>
                    <div class="bg-white p-4 rounded-lg border">
                        <p class="text-sm text-gray-700 mb-2">
                            Instead of guessing, we use <strong>linear algebra</strong> to solve this systematically:
                        </p>
                        <ul class="text-sm text-gray-600 space-y-1 ml-4">
                            <li>• Each element gives us one equation</li>
                            <li>• Each compound gives us one unknown coefficient</li>
                            <li>• We solve the system to get the smallest whole number coefficients</li>
                        </ul>
                    </div>
                </div>
                
                <div class="bg-indigo-50 p-4 rounded-lg border border-indigo-200">
                    <div class="flex items-start">
                        <i class="fas fa-lightbulb text-indigo-600 mt-1 mr-3"></i>
                        <div>
                            <h5 class="font-semibold text-indigo-800 mb-1">Pro Tip</h5>
                            <p class="text-indigo-700 text-sm">
                                In practice, you'd work through this step by step, adjusting coefficients until balanced. 
                                But math gives us the exact answer instantly!
                            </p>
                        </div>
                    </div>
                </div>
            </div>
        `;
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
                        <div class="bg-red-50 p-4 rounded-lg border-l-4 border-red-400">
                            <h5 class="font-bold text-red-700 mb-2">BEFORE (Unbalanced):</h5>
                            <div class="font-mono text-lg text-center bg-white p-3 rounded border">
                                ${result.compounds.map(c => c.formula).join(' + ').replace(/\+ (?=.*[A-Z].*→)/, ' → ')}
                            </div>
                        </div>
                        
                        <div class="flex justify-center">
                            <div class="w-12 h-12 bg-orange-500 text-white rounded-full flex items-center justify-center">
                                <i class="fas fa-arrow-down text-xl"></i>
                            </div>
                        </div>
                        
                        <div class="bg-green-50 p-4 rounded-lg border-l-4 border-green-400">
                            <h5 class="font-bold text-green-700 mb-2">AFTER (Balanced):</h5>
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
                                ${result.compounds.map((compound, i) => {
            const coeff = result.coefficients[i];
            const formatted = coeff === 1 ? compound.formula : `${coeff}${compound.formula}`;
            const color = compound.isReactant ? 'text-red-600' : 'text-green-600';
            return `<span class="${color}">${formatted}</span>`;
        }).join(' + ').replace(/\+ (?=.*text-green)/g, ' → ')}
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

    copyResult() {
        const balancedEquation = document.getElementById('balanced-equation').textContent;
        navigator.clipboard.writeText(balancedEquation).then(() => {
            const copyBtn = document.getElementById('copy-btn');
            const originalText = copyBtn.innerHTML;
            copyBtn.innerHTML = '<i class="fas fa-check mr-1"></i> Copied!';
            copyBtn.classList.remove('bg-gray-600', 'hover:bg-gray-700');
            copyBtn.classList.add('bg-green-600');

            setTimeout(() => {
                copyBtn.innerHTML = originalText;
                copyBtn.classList.remove('bg-green-600');
                copyBtn.classList.add('bg-gray-600', 'hover:bg-gray-700');
            }, 2000);
        });
    }
}

// Initialize the application
document.addEventListener('DOMContentLoaded', () => {
    console.log('DOM loaded, initializing app...');
    new EquationBalancerUI();
});