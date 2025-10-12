class ChemicalEquationPractice {
    constructor() {
        this.currentEquation = null;
        this.correctCount = 0;
        this.incorrectCount = 0;
        this.practiceEquations = [
            // Simple equations
            { reactants: ['H2', 'O2'], products: ['H2O'], coefficients: [2, 1, 2], difficulty: 'easy' },
            { reactants: ['CH4', 'O2'], products: ['CO2', 'H2O'], coefficients: [1, 2, 1, 2], difficulty: 'easy' },
            { reactants: ['N2', 'H2'], products: ['NH3'], coefficients: [1, 3, 2], difficulty: 'easy' },
            { reactants: ['Al', 'O2'], products: ['Al2O3'], coefficients: [4, 3, 2], difficulty: 'easy' },
            { reactants: ['Fe', 'O2'], products: ['Fe2O3'], coefficients: [4, 3, 2], difficulty: 'easy' },
            
            // Medium equations
            { reactants: ['C2H6', 'O2'], products: ['CO2', 'H2O'], coefficients: [2, 7, 4, 6], difficulty: 'medium' },
            { reactants: ['C3H8', 'O2'], products: ['CO2', 'H2O'], coefficients: [1, 5, 3, 4], difficulty: 'medium' },
            { reactants: ['Ca(OH)2', 'HCl'], products: ['CaCl2', 'H2O'], coefficients: [1, 2, 1, 2], difficulty: 'medium' },
            { reactants: ['NaOH', 'H2SO4'], products: ['Na2SO4', 'H2O'], coefficients: [2, 1, 1, 2], difficulty: 'medium' },
            { reactants: ['Mg', 'HCl'], products: ['MgCl2', 'H2'], coefficients: [1, 2, 1, 1], difficulty: 'medium' },
            
            // Hard equations
            { reactants: ['C4H10', 'O2'], products: ['CO2', 'H2O'], coefficients: [2, 13, 8, 10], difficulty: 'hard' },
            { reactants: ['Fe2O3', 'CO'], products: ['Fe', 'CO2'], coefficients: [1, 3, 2, 3], difficulty: 'hard' },
            { reactants: ['KClO3'], products: ['KCl', 'O2'], coefficients: [2, 2, 3], difficulty: 'hard' },
            { reactants: ['Al2(SO4)3', 'Ca(OH)2'], products: ['Al(OH)3', 'CaSO4'], coefficients: [1, 3, 2, 3], difficulty: 'hard' },
            { reactants: ['C6H12O6', 'O2'], products: ['CO2', 'H2O'], coefficients: [1, 6, 6, 6], difficulty: 'hard' }
        ];
        
        this.initializeEventListeners();
        this.updateStats();
    }

    initializeEventListeners() {
        document.getElementById('new-equation-btn').addEventListener('click', () => this.generateNewEquation());
        document.getElementById('check-btn').addEventListener('click', () => this.checkAnswer());
        document.getElementById('hint-btn').addEventListener('click', () => this.showHint());
        document.getElementById('show-solution-btn').addEventListener('click', () => this.showSolution());
        
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
        // Select a random equation
        const randomIndex = Math.floor(Math.random() * this.practiceEquations.length);
        this.currentEquation = { ...this.practiceEquations[randomIndex] };
        
        this.displayEquation();
        this.enableButtons();
        this.clearFeedback();
    }

    displayEquation() {
        const container = document.getElementById('equation-container');
        const { reactants, products } = this.currentEquation;
        
        let html = '<div class="equation-display flex items-center justify-center flex-wrap gap-2">';
        
        // Reactants with input fields
        reactants.forEach((compound, index) => {
            if (index > 0) html += '<span class="text-gray-600 dark:text-gray-400 mx-2">+</span>';
            html += `
                <div class="flex items-center">
                    <input type="text" 
                           class="coefficient-input border-2 border-gray-300 dark:border-gray-600 rounded px-2 py-1 bg-white dark:bg-gray-700 text-gray-800 dark:text-white" 
                           id="coeff-r-${index}" 
                           placeholder="1">
                    <span class="ml-1 text-gray-800 dark:text-white font-mono">${compound}</span>
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
                    <span class="ml-1 text-gray-800 dark:text-white font-mono">${compound}</span>
                </div>
            `;
        });
        
        html += '</div>';
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
            this.showFeedback(true, 'Correct! Well done!');
            this.triggerConfetti();
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

        const hints = [
            "Start by balancing the most complex molecule first.",
            "Count the atoms of each element on both sides.",
            "Try balancing carbon atoms first, then hydrogen, then oxygen.",
            "Remember: coefficients multiply the entire compound.",
            "The sum of atoms for each element must be equal on both sides."
        ];

        const randomHint = hints[Math.floor(Math.random() * hints.length)];
        
        this.showFeedback('hint', `💡 Hint: ${randomHint}`);
    }

    showSolution() {
        if (!this.currentEquation) return;

        const { reactants, products, coefficients } = this.currentEquation;
        let solutionText = '';
        
        // Build solution equation
        reactants.forEach((compound, index) => {
            if (index > 0) solutionText += ' + ';
            solutionText += `${coefficients[index]}${compound}`;
        });
        
        solutionText += ' → ';
        
        products.forEach((compound, index) => {
            if (index > 0) solutionText += ' + ';
            const coeffIndex = reactants.length + index;
            solutionText += `${coefficients[coeffIndex]}${compound}`;
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
    }

    triggerConfetti() {
        // Create confetti container
        const confettiContainer = document.createElement('div');
        confettiContainer.style.position = 'fixed';
        confettiContainer.style.top = '0';
        confettiContainer.style.left = '0';
        confettiContainer.style.width = '100%';
        confettiContainer.style.height = '100%';
        confettiContainer.style.pointerEvents = 'none';
        confettiContainer.style.zIndex = '1000';
        
        // Create multiple confetti pieces
        for (let i = 0; i < 50; i++) {
            const confetti = document.createElement('div');
            confetti.className = 'confetti';
            
            // Random properties for each confetti piece
            const colors = ['#ff6b6b', '#4ecdc4', '#45b7d1', '#96ceb4', '#feca57', '#ff9ff3', '#54a0ff', '#5f27cd', '#00d2d3', '#ff6348', '#2ed573', '#ffa502', '#3742fa', '#70a1ff', '#7bed9f', '#ff7675', '#fd79a8', '#fdcb6e', '#6c5ce7', '#a29bfe'];
            const randomColor = colors[Math.floor(Math.random() * colors.length)];
            
            confetti.style.backgroundColor = randomColor;
            confetti.style.left = Math.random() * 100 + '%';
            confetti.style.animationDelay = Math.random() * 2 + 's';
            confetti.style.animationDuration = (Math.random() * 2 + 2) + 's';
            
            // Random shapes
            if (Math.random() > 0.5) {
                confetti.style.borderRadius = '50%';
            } else {
                confetti.style.transform = 'rotate(45deg)';
            }
            
            confettiContainer.appendChild(confetti);
        }
        
        document.body.appendChild(confettiContainer);
        
        // Remove confetti after animation completes
        setTimeout(() => {
            document.body.removeChild(confettiContainer);
        }, 4000);
    }
}

// Initialize the practice app when the page loads
document.addEventListener('DOMContentLoaded', () => {
    new ChemicalEquationPractice();
});