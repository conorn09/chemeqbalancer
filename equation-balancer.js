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
            
            reactantsInput.addEventListener('focus', function() {
                if (!this.value) {
                    this.placeholder = currentExample.reactants;
                }
            });
            
            productsInput.addEventListener('focus', function() {
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
                title: 'Step 1: Parse the Equation',
                content: this.generateParseStep(originalEquation, result)
            },
            {
                title: 'Step 2: Identify Elements',
                content: this.generateElementsStep(result)
            },
            {
                title: 'Step 3: Set Up System of Equations',
                content: this.generateMatrixStep(result)
            },
            {
                title: 'Step 4: Solve for Coefficients',
                content: this.generateSolutionStep(result)
            },
            {
                title: 'Step 5: Verify Balance',
                content: this.generateVerificationStep(result)
            }
        ];

        steps.forEach((step, index) => {
            const stepDiv = document.createElement('div');
            stepDiv.className = 'mb-4 border border-gray-200 rounded-lg';
            
            stepDiv.innerHTML = `
                <button class="w-full text-left p-4 bg-gray-50 hover:bg-gray-100 transition-colors flex items-center justify-between" 
                        onclick="this.nextElementSibling.classList.toggle('hidden')">
                    <span class="font-medium text-gray-800">${step.title}</span>
                    <i class="fas fa-chevron-down text-gray-500"></i>
                </button>
                <div class="p-4 border-t border-gray-200 ${index === 0 ? '' : 'hidden'}">
                    ${step.content}
                </div>
            `;
            
            container.appendChild(stepDiv);
        });
    }

    generateParseStep(originalEquation, result) {
        const reactants = result.compounds.filter(c => c.isReactant);
        const products = result.compounds.filter(c => !c.isReactant);

        return `
            <p class="mb-3"><strong>Original equation:</strong> <code class="bg-gray-100 px-2 py-1 rounded">${originalEquation}</code></p>
            <p class="mb-3"><strong>Reactants:</strong> ${reactants.map(c => c.formula).join(', ')}</p>
            <p class="mb-3"><strong>Products:</strong> ${products.map(c => c.formula).join(', ')}</p>
            <div class="bg-blue-50 p-3 rounded-lg">
                <p class="text-sm text-blue-800">Each compound is parsed to identify its constituent elements and their quantities.</p>
            </div>
        `;
    }

    generateElementsStep(result) {
        let content = `<p class="mb-3"><strong>Elements found:</strong> ${result.elements.join(', ')}</p>`;
        
        content += '<div class="bg-gray-50 p-3 rounded-lg mb-3"><h4 class="font-medium mb-2">Element breakdown by compound:</h4>';
        result.compounds.forEach(compound => {
            content += `<p class="text-sm"><strong>${compound.formula}:</strong> `;
            const elementList = Object.entries(compound.elements).map(([el, count]) => 
                count === 1 ? el : `${el}×${count}`
            ).join(', ');
            content += elementList + '</p>';
        });
        content += '</div>';

        return content;
    }

    generateMatrixStep(result) {
        let content = '<p class="mb-3">For each element, we create an equation where the sum of atoms on both sides equals zero:</p>';
        
        content += '<div class="bg-gray-50 p-3 rounded-lg mb-3 overflow-x-auto">';
        content += '<table class="w-full text-sm">';
        content += '<thead><tr><th class="text-left p-2">Element</th>';
        
        result.compounds.forEach((compound, i) => {
            content += `<th class="text-center p-2">${compound.formula}</th>`;
        });
        content += '<th class="text-center p-2">=</th></tr></thead><tbody>';

        result.elements.forEach((element, i) => {
            content += `<tr><td class="p-2 font-medium">${element}</td>`;
            result.matrix[i].forEach(coeff => {
                const sign = coeff > 0 ? '+' : '';
                content += `<td class="text-center p-2">${sign}${coeff}</td>`;
            });
            content += '<td class="text-center p-2">0</td></tr>';
        });
        
        content += '</tbody></table></div>';
        
        return content;
    }

    generateSolutionStep(result) {
        let content = '<p class="mb-3">Solving the system of linear equations gives us the coefficients:</p>';
        
        content += '<div class="bg-green-50 p-3 rounded-lg">';
        result.compounds.forEach((compound, i) => {
            content += `<p class="text-sm"><strong>${compound.formula}:</strong> coefficient = ${result.coefficients[i]}</p>`;
        });
        content += '</div>';

        return content;
    }

    generateVerificationStep(result) {
        let content = '<p class="mb-3">Let\'s verify that each element is balanced:</p>';
        
        content += '<div class="bg-gray-50 p-3 rounded-lg">';
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
        content += '</div>';

        content += `<div class="mt-3 p-3 bg-green-50 border border-green-200 rounded-lg">
            <p class="text-green-800 font-medium">✓ Equation is balanced!</p>
            <p class="text-green-700 text-sm mt-1">Final result: <code>${result.balancedEquation}</code></p>
        </div>`;

        return content;
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