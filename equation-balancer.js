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
        const reactantsSearch = document.getElementById('reactants-search');
        const reactantsInput = document.getElementById('reactants-input');

        if (reactantsBtn && reactantsDropdown) {
            reactantsBtn.addEventListener('click', (e) => {
                e.stopPropagation();
                this.toggleDropdown('reactants');
            });

            reactantsInput.addEventListener('focus', () => {
                this.showDropdown('reactants');
            });

            reactantsInput.addEventListener('input', (e) => {
                const query = e.target.value.split('+').pop().trim();
                this.filterDropdown('reactants', query);
            });
        }

        if (reactantsSearch) {
            reactantsSearch.addEventListener('input', (e) => {
                this.filterDropdown('reactants', e.target.value);
            });
        }

        // Products dropdown
        const productsBtn = document.getElementById('products-dropdown-btn');
        const productsDropdown = document.getElementById('products-dropdown');
        const productsSearch = document.getElementById('products-search');
        const productsInput = document.getElementById('products-input');

        if (productsBtn && productsDropdown) {
            productsBtn.addEventListener('click', (e) => {
                e.stopPropagation();
                this.toggleDropdown('products');
            });

            productsInput.addEventListener('focus', () => {
                this.showDropdown('products');
            });

            productsInput.addEventListener('input', (e) => {
                const query = e.target.value.split('+').pop().trim();
                this.filterDropdown('products', query);
            });
        }

        if (productsSearch) {
            productsSearch.addEventListener('input', (e) => {
                this.filterDropdown('products', e.target.value);
            });
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

        const options = listElement.querySelectorAll('.compound-option');
        const categories = listElement.querySelectorAll('.px-3.py-2.text-xs');

        if (!query.trim()) {
            // Show all options
            options.forEach(option => option.style.display = 'block');
            categories.forEach(cat => cat.style.display = 'block');
            return;
        }

        const queryLower = query.toLowerCase();
        let visibleCategories = new Set();

        options.forEach(option => {
            const formula = option.dataset.formula.toLowerCase();
            const name = option.dataset.name.toLowerCase();
            
            if (formula.includes(queryLower) || name.includes(queryLower)) {
                option.style.display = 'block';
                // Find the category for this option
                let categoryElement = option.previousElementSibling;
                while (categoryElement && !categoryElement.classList.contains('text-xs')) {
                    categoryElement = categoryElement.previousElementSibling;
                }
                if (categoryElement) {
                    visibleCategories.add(categoryElement);
                }
            } else {
                option.style.display = 'none';
            }
        });

        // Show/hide categories based on whether they have visible options
        categories.forEach(cat => {
            cat.style.display = visibleCategories.has(cat) ? 'block' : 'none';
        });
    }

    selectCompound(type, formula) {
        const input = document.getElementById(`${type}-input`);
        if (!input) return;

        const currentValue = input.value.trim();
        let newValue;

        if (currentValue === '') {
            newValue = formula;
        } else if (currentValue.endsWith('+')) {
            newValue = currentValue + ' ' + formula;
        } else {
            newValue = currentValue + ' + ' + formula;
        }

        input.value = newValue;
        input.focus();
        this.closeDropdown(type);
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
            
            // Focus search input
            const searchInput = document.getElementById(`${type}-search`);
            if (searchInput) {
                setTimeout(() => searchInput.focus(), 100);
            }
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

        const reactants = reactantsInput.value.trim();
        const products = productsInput.value.trim();

        console.log('Reactants:', reactants);
        console.log('Products:', products);

        if (!reactants || !products) {
            this.showError('Please enter both reactants and products');
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
            this.showResults(result, equation);
        } else {
            this.showError(result.error);
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
        
        // Clear search inputs
        const reactantsSearch = document.getElementById('reactants-search');
        const productsSearch = document.getElementById('products-search');
        if (reactantsSearch) reactantsSearch.value = '';
        if (productsSearch) productsSearch.value = '';
        
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