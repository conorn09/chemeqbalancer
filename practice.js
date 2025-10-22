class ChemicalEquationPractice {
    constructor() {
        this.currentEquation = null;
        this.correctCount = 0;
        this.incorrectCount = 0;
        this.currentDifficulty = 'easy';
        this.completedEquations = {
            easy: new Set(),
            medium: new Set(),
            hard: new Set()
        };
        this.practiceEquations = [
            // Easy equations (12 total)
            { reactants: ['H2', 'O2'], products: ['H2O'], coefficients: [2, 1, 2], difficulty: 'easy' },
            { reactants: ['CH4', 'O2'], products: ['CO2', 'H2O'], coefficients: [1, 2, 1, 2], difficulty: 'easy' },
            { reactants: ['N2', 'H2'], products: ['NH3'], coefficients: [1, 3, 2], difficulty: 'easy' },
            { reactants: ['Al', 'O2'], products: ['Al2O3'], coefficients: [4, 3, 2], difficulty: 'easy' },
            { reactants: ['Fe', 'O2'], products: ['Fe2O3'], coefficients: [4, 3, 2], difficulty: 'easy' },
            { reactants: ['Na', 'Cl2'], products: ['NaCl'], coefficients: [2, 1, 2], difficulty: 'easy' },
            { reactants: ['Mg', 'O2'], products: ['MgO'], coefficients: [2, 1, 2], difficulty: 'easy' },
            { reactants: ['Ca', 'H2O'], products: ['Ca(OH)2', 'H2'], coefficients: [1, 2, 1, 1], difficulty: 'easy' },
            { reactants: ['Li', 'N2'], products: ['Li3N'], coefficients: [6, 1, 2], difficulty: 'easy' },
            { reactants: ['K', 'Br2'], products: ['KBr'], coefficients: [2, 1, 2], difficulty: 'easy' },
            { reactants: ['Zn', 'O2'], products: ['ZnO'], coefficients: [2, 1, 2], difficulty: 'easy' },
            { reactants: ['C', 'O2'], products: ['CO2'], coefficients: [1, 1, 1], difficulty: 'easy' },
            
            // Medium equations (10 total)
            { reactants: ['C2H6', 'O2'], products: ['CO2', 'H2O'], coefficients: [2, 7, 4, 6], difficulty: 'medium' },
            { reactants: ['C3H8', 'O2'], products: ['CO2', 'H2O'], coefficients: [1, 5, 3, 4], difficulty: 'medium' },
            { reactants: ['Ca(OH)2', 'HCl'], products: ['CaCl2', 'H2O'], coefficients: [1, 2, 1, 2], difficulty: 'medium' },
            { reactants: ['NaOH', 'H2SO4'], products: ['Na2SO4', 'H2O'], coefficients: [2, 1, 1, 2], difficulty: 'medium' },
            { reactants: ['Mg', 'HCl'], products: ['MgCl2', 'H2'], coefficients: [1, 2, 1, 1], difficulty: 'medium' },
            { reactants: ['Al', 'HCl'], products: ['AlCl3', 'H2'], coefficients: [2, 6, 2, 3], difficulty: 'medium' },
            { reactants: ['CaCO3', 'HCl'], products: ['CaCl2', 'CO2', 'H2O'], coefficients: [1, 2, 1, 1, 1], difficulty: 'medium' },
            { reactants: ['NH3', 'O2'], products: ['NO', 'H2O'], coefficients: [4, 5, 4, 6], difficulty: 'medium' },
            { reactants: ['Fe', 'H2SO4'], products: ['Fe2(SO4)3', 'H2'], coefficients: [2, 3, 1, 3], difficulty: 'medium' },
            { reactants: ['C2H4', 'O2'], products: ['CO2', 'H2O'], coefficients: [1, 3, 2, 2], difficulty: 'medium' },
            
            // Hard equations (10 total)
            { reactants: ['C4H10', 'O2'], products: ['CO2', 'H2O'], coefficients: [2, 13, 8, 10], difficulty: 'hard' },
            { reactants: ['Fe2O3', 'CO'], products: ['Fe', 'CO2'], coefficients: [1, 3, 2, 3], difficulty: 'hard' },
            { reactants: ['KClO3'], products: ['KCl', 'O2'], coefficients: [2, 2, 3], difficulty: 'hard' },
            { reactants: ['Al2(SO4)3', 'Ca(OH)2'], products: ['Al(OH)3', 'CaSO4'], coefficients: [1, 3, 2, 3], difficulty: 'hard' },
            { reactants: ['C6H12O6', 'O2'], products: ['CO2', 'H2O'], coefficients: [1, 6, 6, 6], difficulty: 'hard' },
            { reactants: ['C8H18', 'O2'], products: ['CO2', 'H2O'], coefficients: [2, 25, 16, 18], difficulty: 'hard' },
            { reactants: ['Ca3(PO4)2', 'H2SO4'], products: ['CaSO4', 'H3PO4'], coefficients: [1, 3, 3, 2], difficulty: 'hard' },
            { reactants: ['Cr2O3', 'Al'], products: ['Al2O3', 'Cr'], coefficients: [1, 2, 1, 2], difficulty: 'hard' },
            { reactants: ['NH4NO3'], products: ['N2O', 'H2O'], coefficients: [1, 1, 2], difficulty: 'hard' },
            { reactants: ['C5H12', 'O2'], products: ['CO2', 'H2O'], coefficients: [1, 8, 5, 6], difficulty: 'hard' }
        ];
        
        this.initializeEventListeners();
        this.updateStats();
    }

    initializeEventListeners() {
        document.getElementById('new-equation-btn').addEventListener('click', () => this.generateNewEquation());
        document.getElementById('check-btn').addEventListener('click', (event) => {
            this.lastClickEvent = event;
            this.checkAnswer();
        });
        document.getElementById('hint-btn').addEventListener('click', () => this.showHint());
        document.getElementById('show-solution-btn').addEventListener('click', () => this.showSolution());
        
        // Difficulty level buttons
        document.getElementById('easy-btn').addEventListener('click', () => this.setDifficulty('easy'));
        document.getElementById('medium-btn').addEventListener('click', () => this.setDifficulty('medium'));
        document.getElementById('hard-btn').addEventListener('click', () => this.setDifficulty('hard'));
        
        // Theme toggle
        const toggleButton = document.getElementById('theme-toggle');
        const themeIcon = document.getElementById('theme-icon');
        const body = document.body;

        // Load saved theme or default to dark mode
        const savedTheme = localStorage.getItem('theme');
        if (savedTheme === 'light') {
            body.classList.remove('dark-mode');
            if (themeIcon) themeIcon.className = 'fas fa-moon text-gray-600';
        } else {
            body.classList.add('dark-mode');
            if (themeIcon) themeIcon.className = 'fas fa-sun text-gray-300';
            if (!savedTheme) localStorage.setItem('theme', 'dark');
        }

        if (toggleButton) {
            toggleButton.addEventListener('click', () => {
                body.classList.toggle('dark-mode');
                const isDark = body.classList.contains('dark-mode');
                if (themeIcon) {
                    themeIcon.className = isDark ? 'fas fa-sun text-gray-300' : 'fas fa-moon text-gray-600';
                }
                localStorage.setItem('theme', isDark ? 'dark' : 'light');
            });
        }
    }

    generateNewEquation() {
        // Filter equations by current difficulty
        const filteredEquations = this.practiceEquations.filter(eq => eq.difficulty === this.currentDifficulty);
        
        // Get equations that haven't been completed yet
        const availableEquations = filteredEquations.filter((eq, index) => {
            const equationId = this.getEquationId(eq);
            return !this.completedEquations[this.currentDifficulty].has(equationId);
        });
        
        // If all equations are completed, reset the completed list and show congratulations
        if (availableEquations.length === 0) {
            this.completedEquations[this.currentDifficulty].clear();
            this.showFeedback('solution', `🎉 Congratulations! You've completed all ${this.currentDifficulty} level equations! Starting over with fresh questions.`);
            // Use all equations again
            const resetEquations = filteredEquations;
            const randomIndex = Math.floor(Math.random() * resetEquations.length);
            this.currentEquation = { ...resetEquations[randomIndex] };
        } else {
            // Select a random equation from available (not completed) equations
            const randomIndex = Math.floor(Math.random() * availableEquations.length);
            this.currentEquation = { ...availableEquations[randomIndex] };
        }
        
        this.displayEquation();
        this.enableButtons();
        this.clearFeedback();
    }

    getEquationId(equation) {
        // Create a unique identifier for each equation
        return equation.reactants.join('+') + '->' + equation.products.join('+');
    }

    setDifficulty(difficulty) {
        this.currentDifficulty = difficulty;
        
        // Update button styles
        document.querySelectorAll('.difficulty-btn').forEach(btn => {
            btn.classList.remove('bg-emerald-500', 'bg-amber-500', 'bg-red-500', 'text-white');
            btn.classList.add('text-gray-600', 'dark:text-gray-300', 'hover:bg-gray-200', 'dark:hover:bg-gray-600');
        });
        
        // Highlight selected difficulty
        const selectedBtn = document.getElementById(`${difficulty}-btn`);
        selectedBtn.classList.remove('text-gray-600', 'dark:text-gray-300', 'hover:bg-gray-200', 'dark:hover:bg-gray-600');
        
        if (difficulty === 'easy') {
            selectedBtn.classList.add('bg-emerald-500', 'text-white');
        } else if (difficulty === 'medium') {
            selectedBtn.classList.add('bg-amber-500', 'text-white');
        } else if (difficulty === 'hard') {
            selectedBtn.classList.add('bg-red-500', 'text-white');
        }
        
        // Update button text to show progress
        this.updateDifficultyProgress();
        
        // Generate new equation with selected difficulty
        this.generateNewEquation();
    }

    updateDifficultyProgress() {
        const difficulties = ['easy', 'medium', 'hard'];
        
        difficulties.forEach(diff => {
            const btn = document.getElementById(`${diff}-btn`);
            const totalQuestions = this.practiceEquations.filter(eq => eq.difficulty === diff).length;
            const completedCount = this.completedEquations[diff].size;
            
            // Update button text to show progress
            const icon = diff === 'easy' ? 'fas fa-seedling' : diff === 'medium' ? 'fas fa-leaf' : 'fas fa-fire';
            const diffName = diff.charAt(0).toUpperCase() + diff.slice(1);
            
            if (completedCount > 0) {
                btn.innerHTML = `<i class="${icon} mr-1"></i>${diffName} (${completedCount}/${totalQuestions})`;
            } else {
                btn.innerHTML = `<i class="${icon} mr-1"></i>${diffName}`;
            }
        });
    }

    displayEquation() {
        const container = document.getElementById('equation-container');
        const { reactants, products } = this.currentEquation;
        
        // Add difficulty indicator
        const difficultyColors = {
            easy: 'bg-emerald-100 text-emerald-800 dark:bg-emerald-900 dark:text-emerald-200',
            medium: 'bg-amber-100 text-amber-800 dark:bg-amber-900 dark:text-amber-200',
            hard: 'bg-red-100 text-red-800 dark:bg-red-900 dark:text-red-200'
        };
        
        let html = `
            <div class="flex items-center justify-center mb-3">
                <span class="px-3 py-1 rounded-full text-xs font-medium ${difficultyColors[this.currentDifficulty]}">
                    ${this.currentDifficulty.charAt(0).toUpperCase() + this.currentDifficulty.slice(1)} Level
                </span>
            </div>
            <div class="equation-display flex items-center justify-center flex-wrap gap-2">
        `;
        
        // Reactants with input fields
        reactants.forEach((compound, index) => {
            if (index > 0) html += '<span class="text-gray-600 dark:text-gray-400 mx-2">+</span>';
            html += `
                <div class="flex items-center">
                    <input type="text" 
                           class="coefficient-input border-2 border-gray-300 dark:border-gray-600 rounded px-2 py-1 bg-white dark:bg-gray-700 text-gray-800 dark:text-white" 
                           id="coeff-r-${index}" 
                           placeholder="1">
                    <span class="ml-1 text-gray-800 dark:text-white font-mono">${this.formatChemicalFormula(compound)}</span>
                </div>
            `;
        });
        
        html += '<span class="text-gray-600 dark:text-gray-400 mx-4 text-2xl">→</span>';
        
        // Products with input fields
        products.forEach((compound, index) => {
            if (index > 0) html += '<span class="text-gray-600 dark:text-gray-400 mx-2">+</span>';
            html += `
                <div class="flex items-center">
                    <input type="text" 
                           class="coefficient-input border-2 border-gray-300 dark:border-gray-600 rounded px-2 py-1 bg-white dark:bg-gray-700 text-gray-800 dark:text-white" 
                           id="coeff-p-${index}" 
                           placeholder="1">
                    <span class="ml-1 text-gray-800 dark:text-white font-mono">${this.formatChemicalFormula(compound)}</span>
                </div>
            `;
        });
        
        html += '</div></div>';
        container.innerHTML = html;

        // Add enter key listener and input validation to inputs
        container.querySelectorAll('input').forEach(input => {
            input.addEventListener('keypress', (e) => {
                if (e.key === 'Enter') {
                    this.checkAnswer();
                }
                // Only allow numbers
                if (!/[0-9]/.test(e.key) && e.key !== 'Backspace' && e.key !== 'Delete' && e.key !== 'Enter') {
                    e.preventDefault();
                }
            });
            
            // Validate input on blur
            input.addEventListener('input', (e) => {
                const value = e.target.value;
                if (value && (isNaN(value) || parseInt(value) <= 0)) {
                    e.target.value = '';
                }
            });
        });
    }

    enableButtons() {
        document.getElementById('check-btn').disabled = false;
        document.getElementById('hint-btn').disabled = false;
        document.getElementById('show-solution-btn').disabled = false;
    }

    checkAnswer() {
        if (!this.currentEquation) return;

        const userCoefficients = this.getUserCoefficients();
        const correctCoefficients = this.currentEquation.coefficients;
        
        // Check if any coefficients are invalid (0 or negative)
        if (userCoefficients.some(coeff => coeff <= 0)) {
            this.showFeedback(false, 'Please enter only positive numbers for coefficients.');
            return;
        }

        // Check if coefficients are correct (allowing for multiples)
        const isCorrect = this.areCoefficientsEquivalent(userCoefficients, correctCoefficients);
        
        if (isCorrect) {
            this.correctCount++;
            
            // Mark this equation as completed
            const equationId = this.getEquationId(this.currentEquation);
            this.completedEquations[this.currentDifficulty].add(equationId);
            
            this.showFeedback(true, 'Correct! Well done!');
            // Store the click event for confetti explosion
            this.triggerConfetti(this.lastClickEvent);
        } else {
            this.incorrectCount++;
            this.showFeedback(false, 'Incorrect. Try again or use the hint button for help.');
        }
        
        this.updateStats();
    }

    getUserCoefficients() {
        const coefficients = [];
        const { reactants, products } = this.currentEquation;
        
        // Get reactant coefficients
        for (let i = 0; i < reactants.length; i++) {
            const input = document.getElementById(`coeff-r-${i}`);
            const value = parseInt(input.value) || 1; // Default to 1 if empty
            coefficients.push(value);
        }
        
        // Get product coefficients
        for (let i = 0; i < products.length; i++) {
            const input = document.getElementById(`coeff-p-${i}`);
            const value = parseInt(input.value) || 1; // Default to 1 if empty
            coefficients.push(value);
        }
        
        return coefficients;
    }

    areCoefficientsEquivalent(userCoeffs, correctCoeffs) {
        // Find the GCD to normalize both sets of coefficients
        const gcd = (a, b) => b === 0 ? a : gcd(b, a % b);
        const arrayGcd = arr => arr.reduce(gcd);
        
        const userGcd = arrayGcd(userCoeffs.filter(c => c > 0));
        const correctGcd = arrayGcd(correctCoeffs);
        
        const normalizedUser = userCoeffs.map(c => c / userGcd);
        const normalizedCorrect = correctCoeffs.map(c => c / correctGcd);
        
        return normalizedUser.every((coeff, index) => coeff === normalizedCorrect[index]);
    }

    showHint() {
        if (!this.currentEquation) return;

        const hintsByDifficulty = {
            easy: [
                "Start with the element that appears in the fewest compounds.",
                "Count atoms on each side - they must be equal.",
                "Try small numbers first: 1, 2, 3, 4.",
                "Balance one element at a time.",
                "Remember: coefficients go in front of the entire compound."
            ],
            medium: [
                "Start by balancing the most complex molecule first.",
                "Try balancing carbon atoms first, then hydrogen, then oxygen.",
                "Look for patterns - some coefficients might be related.",
                "If you get fractions, multiply everything to get whole numbers.",
                "Count carefully - polyatomic ions can be tricky."
            ],
            hard: [
                "Use the algebraic method for complex equations.",
                "Balance metals first, then non-metals, then H and O last.",
                "Consider using the least common multiple approach.",
                "Some equations may require coefficients larger than 10.",
                "Break down complex polyatomic ions step by step."
            ]
        };

        const hints = hintsByDifficulty[this.currentDifficulty];
        const randomHint = hints[Math.floor(Math.random() * hints.length)];
        
        this.showFeedback('hint', `💡 ${this.currentDifficulty.toUpperCase()} Hint: ${randomHint}`);
    }

    showSolution() {
        if (!this.currentEquation) return;

        const { reactants, products, coefficients } = this.currentEquation;
        let solutionText = '';
        
        // Build solution equation with proper subscripts
        reactants.forEach((compound, index) => {
            if (index > 0) solutionText += ' + ';
            solutionText += `${coefficients[index]}${this.formatChemicalFormula(compound)}`;
        });
        
        solutionText += ' → ';
        
        products.forEach((compound, index) => {
            if (index > 0) solutionText += ' + ';
            const coeffIndex = reactants.length + index;
            solutionText += `${coefficients[coeffIndex]}${this.formatChemicalFormula(compound)}`;
        });

        // Fill in the correct coefficients
        reactants.forEach((_, index) => {
            const input = document.getElementById(`coeff-r-${index}`);
            input.value = coefficients[index];
        });
        
        products.forEach((_, index) => {
            const input = document.getElementById(`coeff-p-${index}`);
            const coeffIndex = reactants.length + index;
            input.value = coefficients[coeffIndex];
        });

        this.showFeedback('solution', `✅ Solution: ${solutionText}`);
    }

    showFeedback(type, message) {
        const container = document.getElementById('feedback-container');
        container.classList.remove('hidden');
        
        let bgColor, textColor, icon;
        
        switch (type) {
            case true: // correct
                bgColor = 'bg-green-50 dark:bg-green-900 border-green-200 dark:border-green-700';
                textColor = 'text-green-800 dark:text-green-200';
                icon = 'fas fa-check-circle text-green-600 dark:text-green-400';
                container.classList.add('feedback-correct');
                break;
            case false: // incorrect
                bgColor = 'bg-red-50 dark:bg-red-900 border-red-200 dark:border-red-700';
                textColor = 'text-red-800 dark:text-red-200';
                icon = 'fas fa-times-circle text-red-600 dark:text-red-400';
                container.classList.add('feedback-incorrect');
                break;
            case 'hint':
                bgColor = 'bg-amber-50 dark:bg-amber-900 border-amber-200 dark:border-amber-700';
                textColor = 'text-amber-800 dark:text-amber-200';
                icon = 'fas fa-lightbulb text-amber-600 dark:text-amber-400';
                break;
            case 'solution':
                bgColor = 'bg-blue-50 dark:bg-blue-900 border-blue-200 dark:border-blue-700';
                textColor = 'text-blue-800 dark:text-blue-200';
                icon = 'fas fa-eye text-blue-600 dark:text-blue-400';
                break;
        }
        
        container.innerHTML = `
            <div class="p-4 rounded-lg border-2 ${bgColor} ${textColor} mb-6">
                <div class="flex items-center">
                    <i class="${icon} mr-3"></i>
                    <span class="font-medium">${message}</span>
                </div>
            </div>
        `;

        // Remove animation classes after animation completes
        setTimeout(() => {
            container.classList.remove('feedback-correct', 'feedback-incorrect');
        }, 600);
    }

    clearFeedback() {
        const container = document.getElementById('feedback-container');
        container.classList.add('hidden');
        container.innerHTML = '';
    }

    updateStats() {
        document.getElementById('correct-count').textContent = this.correctCount;
        document.getElementById('incorrect-count').textContent = this.incorrectCount;
        
        const total = this.correctCount + this.incorrectCount;
        const accuracy = total > 0 ? Math.round((this.correctCount / total) * 100) : 0;
        document.getElementById('accuracy').textContent = `${accuracy}%`;
        
        // Update difficulty progress
        this.updateDifficultyProgress();
    }

    formatChemicalFormula(formula) {
        // Convert numbers in chemical formulas to subscripts
        return formula.replace(/(\d+)/g, '<sub>$1</sub>');
    }

    triggerConfetti(clickEvent) {
        let explosionX, explosionY;
        
        if (clickEvent && clickEvent.clientX && clickEvent.clientY) {
            // Use exact click position
            explosionX = clickEvent.clientX;
            explosionY = clickEvent.clientY;
        } else {
            // Fallback to button center
            const checkButton = document.getElementById('check-btn');
            const buttonRect = checkButton.getBoundingClientRect();
            explosionX = buttonRect.left + buttonRect.width / 2;
            explosionY = buttonRect.top + buttonRect.height / 2;
        }
        
        const colors = ['#ff6b6b', '#4ecdc4', '#45b7d1', '#96ceb4', '#feca57', '#ff9ff3', '#54a0ff', '#5f27cd', '#00d2d3', '#ff6348', '#2ed573', '#ffa502', '#3742fa', '#70a1ff', '#7bed9f', '#ff7675', '#fd79a8', '#fdcb6e', '#6c5ce7', '#a29bfe'];
        
        // Create initial explosion (thinner)
        for (let wave = 0; wave < 3; wave++) {
            setTimeout(() => {
                for (let i = 0; i < 30; i++) {
                    const confetti = document.createElement('div');
                    confetti.className = 'confetti';
                    
                    // Random color
                    const randomColor = colors[Math.floor(Math.random() * colors.length)];
                    confetti.style.backgroundColor = randomColor;
                    
                    // Random shape and size
                    if (Math.random() > 0.5) {
                        confetti.style.borderRadius = '50%';
                    }
                    
                    // Vary size more dramatically for better visibility
                    const size = 8 + Math.random() * 8; // Larger pieces (8-16px)
                    confetti.style.width = size + 'px';
                    confetti.style.height = size + 'px';
                    
                    // Start at explosion position with slight random offset
                    const offsetX = (Math.random() - 0.5) * 30;
                    const offsetY = (Math.random() - 0.5) * 30;
                    confetti.style.left = (explosionX + offsetX) + 'px';
                    confetti.style.top = (explosionY + offsetY) + 'px';
                    
                    // Random explosion direction and speed - make it MASSIVE
                    const angle = Math.random() * Math.PI * 2; // Full 360 degree spread
                    const velocity = 25 + Math.random() * 35; // Extremely high initial velocity (25-60) for maximum coverage
                    const velocityX = Math.cos(angle) * velocity;
                    const velocityY = Math.sin(angle) * velocity;
                    
                    document.body.appendChild(confetti);
                    
                    // Animate the confetti piece
                    this.animateConfettiPiece(confetti, velocityX, velocityY, explosionX + offsetX, explosionY + offsetY);
                }
            }, wave * 50); // Stagger waves by 50ms
        }
        
        // Create dense falling confetti shower from the top
        this.createFallingConfettiShower(colors);
    }

    animateConfettiPiece(confetti, velocityX, velocityY, startX, startY) {
        let x = startX;
        let y = startY;
        let vx = velocityX;
        let vy = velocityY;
        let rotation = 0;
        let opacity = 1;
        
        const gravity = 0.25; // Much slower falling
        const friction = 0.995; // Less air resistance for longer flight
        const rotationSpeed = (Math.random() - 0.5) * 20;
        
        const animate = () => {
            // Apply physics
            vy += gravity; // Gravity
            vx *= friction; // Air resistance
            vy *= friction;
            
            // Update position
            x += vx;
            y += vy;
            rotation += rotationSpeed;
            
            // Fade out over time (very slow fade for maximum visibility)
            opacity -= 0.002;
            
            // Update confetti position and rotation
            confetti.style.left = x + 'px';
            confetti.style.top = y + 'px';
            confetti.style.transform = `rotate(${rotation}deg)`;
            confetti.style.opacity = opacity;
            
            // Continue animation if confetti is still visible and on screen
            if (opacity > 0 && y < window.innerHeight + 50) {
                requestAnimationFrame(animate);
            } else {
                // Remove confetti when it's done
                if (confetti.parentNode) {
                    confetti.parentNode.removeChild(confetti);
                }
            }
        };
        
        requestAnimationFrame(animate);
    }

    createFallingConfettiShower(colors) {
        // Create dense falling confetti from the top of the screen
        const screenWidth = window.innerWidth;
        
        // Generate confetti continuously for 3 seconds
        let confettiInterval = setInterval(() => {
            // Create a burst of confetti across the top of the screen
            for (let i = 0; i < 15; i++) {
                const confetti = document.createElement('div');
                confetti.className = 'confetti';
                
                // Random color
                const randomColor = colors[Math.floor(Math.random() * colors.length)];
                confetti.style.backgroundColor = randomColor;
                
                // Random shape and size
                if (Math.random() > 0.5) {
                    confetti.style.borderRadius = '50%';
                }
                
                const size = 6 + Math.random() * 8;
                confetti.style.width = size + 'px';
                confetti.style.height = size + 'px';
                
                // Start from random position across the top of the screen
                const startX = Math.random() * screenWidth;
                const startY = -20; // Start above the screen
                
                confetti.style.left = startX + 'px';
                confetti.style.top = startY + 'px';
                
                // Falling motion with slight horizontal drift
                const velocityX = (Math.random() - 0.5) * 4; // Slight horizontal drift
                const velocityY = 2 + Math.random() * 3; // Downward velocity
                
                document.body.appendChild(confetti);
                
                // Animate the falling confetti
                this.animateFallingConfetti(confetti, velocityX, velocityY, startX, startY);
            }
        }, 100); // Create new confetti every 100ms
        
        // Stop creating new confetti after 3 seconds
        setTimeout(() => {
            clearInterval(confettiInterval);
        }, 3000);
    }

    animateFallingConfetti(confetti, velocityX, velocityY, startX, startY) {
        let x = startX;
        let y = startY;
        let vx = velocityX;
        let vy = velocityY;
        let rotation = 0;
        let opacity = 1;
        
        const gravity = 0.1; // Light gravity for gentle falling
        const rotationSpeed = (Math.random() - 0.5) * 10;
        
        const animate = () => {
            // Apply gentle physics
            vy += gravity;
            
            // Update position
            x += vx;
            y += vy;
            rotation += rotationSpeed;
            
            // Fade out slowly
            opacity -= 0.002;
            
            // Update confetti position and rotation
            confetti.style.left = x + 'px';
            confetti.style.top = y + 'px';
            confetti.style.transform = `rotate(${rotation}deg)`;
            confetti.style.opacity = opacity;
            
            // Continue animation if confetti is still visible and on screen
            if (opacity > 0 && y < window.innerHeight + 100) {
                requestAnimationFrame(animate);
            } else {
                // Remove confetti when it's done
                if (confetti.parentNode) {
                    confetti.parentNode.removeChild(confetti);
                }
            }
        };
        
        requestAnimationFrame(animate);
    }
}

// Initialize the practice app when the page loads
document.addEventListener('DOMContentLoaded', () => {
    new ChemicalEquationPractice();
});