class AtomicTheoryPractice {
    constructor() {
        this.currentQuestion = null;
        this.correctCount = 0;
        this.incorrectCount = 0;
        this.hintTimer = null;
        this.lastActivity = Date.now();
        
        // Element data with atomic numbers and common isotopes
        this.elements = [
            { symbol: 'H', name: 'Hydrogen', atomicNumber: 1, commonMasses: [1, 2, 3] },
            { symbol: 'He', name: 'Helium', atomicNumber: 2, commonMasses: [3, 4] },
            { symbol: 'Li', name: 'Lithium', atomicNumber: 3, commonMasses: [6, 7] },
            { symbol: 'Be', name: 'Beryllium', atomicNumber: 4, commonMasses: [9] },
            { symbol: 'B', name: 'Boron', atomicNumber: 5, commonMasses: [10, 11] },
            { symbol: 'C', name: 'Carbon', atomicNumber: 6, commonMasses: [12, 13, 14] },
            { symbol: 'N', name: 'Nitrogen', atomicNumber: 7, commonMasses: [14, 15] },
            { symbol: 'O', name: 'Oxygen', atomicNumber: 8, commonMasses: [16, 17, 18] },
            { symbol: 'F', name: 'Fluorine', atomicNumber: 9, commonMasses: [19] },
            { symbol: 'Ne', name: 'Neon', atomicNumber: 10, commonMasses: [20, 21, 22] },
            { symbol: 'Na', name: 'Sodium', atomicNumber: 11, commonMasses: [23] },
            { symbol: 'Mg', name: 'Magnesium', atomicNumber: 12, commonMasses: [24, 25, 26] },
            { symbol: 'Al', name: 'Aluminum', atomicNumber: 13, commonMasses: [27] },
            { symbol: 'Si', name: 'Silicon', atomicNumber: 14, commonMasses: [28, 29, 30] },
            { symbol: 'P', name: 'Phosphorus', atomicNumber: 15, commonMasses: [31] },
            { symbol: 'S', name: 'Sulfur', atomicNumber: 16, commonMasses: [32, 33, 34, 36] },
            { symbol: 'Cl', name: 'Chlorine', atomicNumber: 17, commonMasses: [35, 37] },
            { symbol: 'Ar', name: 'Argon', atomicNumber: 18, commonMasses: [36, 38, 40] },
            { symbol: 'K', name: 'Potassium', atomicNumber: 19, commonMasses: [39, 40, 41] },
            { symbol: 'Ca', name: 'Calcium', atomicNumber: 20, commonMasses: [40, 42, 43, 44, 46, 48] }
        ];
        
        this.initializeEventListeners();
        this.updateStats();
    }

    initializeEventListeners() {
        document.getElementById('new-question-btn').addEventListener('click', () => this.generateNewQuestion());
        document.getElementById('check-answer-btn').addEventListener('click', () => this.checkAnswer());
        document.getElementById('hint-btn').addEventListener('click', () => this.showHint());
        document.getElementById('show-answer-btn').addEventListener('click', () => this.showAnswer());
        
        // Add enter key listeners to inputs
        ['protons-input', 'neutrons-input', 'electrons-input'].forEach(id => {
            const input = document.getElementById(id);
            input.addEventListener('keypress', (e) => {
                if (e.key === 'Enter') {
                    this.checkAnswer();
                }
                // Only allow numbers
                if (!/[0-9]/.test(e.key) && e.key !== 'Backspace' && e.key !== 'Delete' && e.key !== 'Enter') {
                    e.preventDefault();
                }
            });
            
            // Validate input and track activity
            input.addEventListener('input', (e) => {
                const value = e.target.value;
                if (value && (isNaN(value) || parseInt(value) < 0)) {
                    e.target.value = '';
                }
                this.trackActivity();
            });
        });
        
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

    generateNewQuestion() {
        // Select random element
        const element = this.elements[Math.floor(Math.random() * this.elements.length)];
        
        // Select random mass number from common isotopes
        const massNumber = element.commonMasses[Math.floor(Math.random() * element.commonMasses.length)];
        
        // Randomly decide if it's an ion (30% chance)
        const isIon = Math.random() < 0.3;
        let charge = 0;
        
        if (isIon) {
            // Generate realistic ion charges based on element
            const commonCharges = this.getCommonCharges(element.atomicNumber);
            charge = commonCharges[Math.floor(Math.random() * commonCharges.length)];
        }
        
        this.currentQuestion = {
            element: element,
            atomicNumber: element.atomicNumber,
            massNumber: massNumber,
            charge: charge,
            protons: element.atomicNumber,
            neutrons: massNumber - element.atomicNumber,
            electrons: element.atomicNumber - charge
        };
        
        this.displayQuestion();
        this.showAnswerInputs();
        this.clearInputs();
        this.clearFeedback();
        this.startHintTimer();
    }

    trackActivity() {
        this.lastActivity = Date.now();
        this.stopHintAnimation();
        this.startHintTimer();
    }

    startHintTimer() {
        // Clear any existing timer
        if (this.hintTimer) {
            clearTimeout(this.hintTimer);
        }
        
        // Start new timer for 10 seconds
        this.hintTimer = setTimeout(() => {
            this.animateHintButton();
        }, 10000);
    }

    stopHintAnimation() {
        const hintBtn = document.getElementById('hint-btn');
        if (hintBtn) {
            hintBtn.classList.remove('hint-bounce', 'hint-wiggle');
        }
    }

    animateHintButton() {
        const hintBtn = document.getElementById('hint-btn');
        if (hintBtn && this.currentQuestion) {
            // First do a wiggle, then start bouncing
            hintBtn.classList.add('hint-wiggle');
            
            setTimeout(() => {
                hintBtn.classList.remove('hint-wiggle');
                hintBtn.classList.add('hint-bounce');
                
                // Stop bouncing after 5 seconds
                setTimeout(() => {
                    hintBtn.classList.remove('hint-bounce');
                }, 5000);
            }, 1500); // After wiggle animation completes
        }
    }

    getCommonCharges(atomicNumber) {
        // Return common ion charges based on element position
        const chargeMap = {
            1: [-1, 1], // H
            3: [1], // Li
            4: [2], // Be
            5: [3], // B
            6: [-4, 4], // C
            7: [-3, 3, 5], // N
            8: [-2], // O
            9: [-1], // F
            11: [1], // Na
            12: [2], // Mg
            13: [3], // Al
            14: [-4, 4], // Si
            15: [-3, 3, 5], // P
            16: [-2, 4, 6], // S
            17: [-1, 1, 3, 5, 7], // Cl
            19: [1], // K
            20: [2] // Ca
        };
        
        return chargeMap[atomicNumber] || [1, -1, 2, -2];
    }

    displayQuestion() {
        const container = document.getElementById('question-container');
        const { element, atomicNumber, massNumber, charge } = this.currentQuestion;
        
        const chargeDisplay = charge === 0 ? '' : 
            charge > 0 ? `<sup>${charge > 1 ? charge : ''}+</sup>` : 
            `<sup>${Math.abs(charge) > 1 ? Math.abs(charge) : ''}-</sup>`;
        
        container.innerHTML = `
            <div class="text-center">
                <div class="mb-4">
                    <h3 class="text-2xl font-bold text-gray-800 dark:text-white mb-4">
                        Determine the number of protons, neutrons, and electrons
                    </h3>
                </div>
                
                <div class="bg-white dark:bg-gray-600 rounded-lg p-6 border-2 border-blue-300 dark:border-blue-600 inline-block">
                    <div class="text-center">
                        <div class="text-4xl font-bold text-blue-600 dark:text-blue-300 mb-2">
                            <sup>${massNumber}</sup>${element.symbol}${chargeDisplay}
                        </div>
                        <div class="text-sm text-gray-600 dark:text-gray-300 space-y-1">
                            <div><strong>Element:</strong> ${element.name}</div>
                            <div><strong>Atomic Number:</strong> ${atomicNumber}</div>
                            <div><strong>Mass Number:</strong> ${massNumber}</div>
                            ${charge !== 0 ? `<div><strong>Charge:</strong> ${charge > 0 ? '+' : ''}${charge}</div>` : ''}
                        </div>
                    </div>
                </div>
                

            </div>
        `;
    }

    showAnswerInputs() {
        document.getElementById('answer-container').classList.remove('hidden');
    }

    clearInputs() {
        document.getElementById('protons-input').value = '';
        document.getElementById('neutrons-input').value = '';
        document.getElementById('electrons-input').value = '';
    }

    checkAnswer() {
        if (!this.currentQuestion) return;

        // Stop hint animation when user checks answer
        this.stopHintAnimation();
        if (this.hintTimer) {
            clearTimeout(this.hintTimer);
        }

        const userProtons = parseInt(document.getElementById('protons-input').value) || 0;
        const userNeutrons = parseInt(document.getElementById('neutrons-input').value) || 0;
        const userElectrons = parseInt(document.getElementById('electrons-input').value) || 0;

        const { protons, neutrons, electrons, element, massNumber } = this.currentQuestion;

        const protonsCorrect = userProtons === protons;
        const neutronsCorrect = userNeutrons === neutrons;
        const electronsCorrect = userElectrons === electrons;

        const allCorrect = protonsCorrect && neutronsCorrect && electronsCorrect;

        if (allCorrect) {
            this.correctCount++;
            this.showFeedback(true, `✅ Correct! That's ${element.name}-${massNumber}.`);
        } else {
            this.incorrectCount++;
            let errorMsg = '❌ Not quite right. ';
            const errors = [];
            
            if (!protonsCorrect) errors.push(`Protons should be ${protons}`);
            if (!neutronsCorrect) errors.push(`Neutrons should be ${neutrons}`);
            if (!electronsCorrect) errors.push(`Electrons should be ${electrons}`);
            
            errorMsg += errors.join(', ') + '.';
            this.showFeedback(false, errorMsg);
        }

        this.updateStats();
    }

    showHint() {
        if (!this.currentQuestion) return;

        // Stop hint animation when clicked
        this.stopHintAnimation();
        if (this.hintTimer) {
            clearTimeout(this.hintTimer);
        }

        const hints = [
            "The atomic number tells you the number of protons.",
            "Mass number = protons + neutrons, so neutrons = mass number - protons.",
            "In a neutral atom, electrons = protons. For ions, electrons = protons - charge.",
            "Positive ions have fewer electrons than protons.",
            "Negative ions have more electrons than protons.",
            "Remember: Atomic number = protons, Mass number = protons + neutrons.",
            "For ions: Electrons = protons - charge (positive charge means fewer electrons)."
        ];

        const randomHint = hints[Math.floor(Math.random() * hints.length)];
        this.showFeedback('hint', `💡 Hint: ${randomHint}`);
    }

    showAnswer() {
        if (!this.currentQuestion) return;

        const { protons, neutrons, electrons, element, massNumber } = this.currentQuestion;

        // Fill in the correct answers
        document.getElementById('protons-input').value = protons;
        document.getElementById('neutrons-input').value = neutrons;
        document.getElementById('electrons-input').value = electrons;

        this.showFeedback('solution', `✅ Answer: ${protons} protons, ${neutrons} neutrons, ${electrons} electrons for ${element.name}-${massNumber}.`);
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
}

// Initialize the atomic practice app when the page loads
document.addEventListener('DOMContentLoaded', () => {
    new AtomicTheoryPractice();
});