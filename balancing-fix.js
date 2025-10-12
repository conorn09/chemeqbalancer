// Fixed simulateBalancingProcess function
function simulateBalancingProcess(result, balancingOrder) {
    const steps = [];

    // Start with all coefficients as 1 (unbalanced)
    let currentCoeffs = new Array(result.compounds.length).fill(1);

    balancingOrder.forEach((element, index) => {
        // BEFORE balancing this element - show current state
        const beforeCoeffs = [...currentCoeffs];

        // Calculate what the coefficients should be AFTER balancing this element
        const afterCoeffs = [...currentCoeffs];

        // Only update coefficients for compounds that contain the current element
        result.compounds.forEach((compound, i) => {
            if (compound.elements[element]) {
                afterCoeffs[i] = result.coefficients[i];
            }
        });

        // Calculate atom counts BEFORE balancing this element
        let beforeLeftCount = 0;
        let beforeRightCount = 0;
        let beforeLeftCompounds = [];
        let beforeRightCompounds = [];

        result.compounds.forEach((compound, i) => {
            const elementCount = compound.elements[element] || 0;
            const coeff = beforeCoeffs[i];
            const totalCount = elementCount * coeff;

            const compoundDisplay = coeff === 1 ? compound.formula : `${coeff}${compound.formula}`;
            const elementBreakdown = elementCount > 0 ?
                `Contains ${elementCount} ${element} atom${elementCount > 1 ? 's' : ''} × ${coeff} = ${totalCount} total` :
                `No ${element} atoms`;

            if (compound.isReactant) {
                beforeLeftCount += totalCount;
                if (elementCount > 0) {
                    beforeLeftCompounds.push({
                        display: compoundDisplay,
                        elementBreakdown: elementBreakdown
                    });
                }
            } else {
                beforeRightCount += totalCount;
                if (elementCount > 0) {
                    beforeRightCompounds.push({
                        display: compoundDisplay,
                        elementBreakdown: elementBreakdown
                    });
                }
            }
        });

        // Calculate atom counts AFTER balancing this element
        let afterLeftCount = 0;
        let afterRightCount = 0;
        let afterLeftCompounds = [];
        let afterRightCompounds = [];

        result.compounds.forEach((compound, i) => {
            const elementCount = compound.elements[element] || 0;
            const coeff = afterCoeffs[i];
            const totalCount = elementCount * coeff;

            const compoundDisplay = coeff === 1 ? compound.formula : `${coeff}${compound.formula}`;
            const elementBreakdown = elementCount > 0 ?
                `Contains ${elementCount} ${element} atom${elementCount > 1 ? 's' : ''} × ${coeff} = ${totalCount} total` :
                `No ${element} atoms`;

            if (compound.isReactant) {
                afterLeftCount += totalCount;
                if (elementCount > 0) {
                    afterLeftCompounds.push({
                        display: compoundDisplay,
                        elementBreakdown: elementBreakdown
                    });
                }
            } else {
                afterRightCount += totalCount;
                if (elementCount > 0) {
                    afterRightCompounds.push({
                        display: compoundDisplay,
                        elementBreakdown: elementBreakdown
                    });
                }
            }
        });

        // Generate equation states
        const beforeEquation = result.compounds.map((compound, i) => {
            const coeff = beforeCoeffs[i];
            return coeff === 1 ? compound.formula : `${coeff}${compound.formula}`;
        }).join(' + ').replace(/\+ (?=.*[A-Z].*→)/, ' → ');

        const afterEquation = result.compounds.map((compound, i) => {
            const coeff = afterCoeffs[i];
            return coeff === 1 ? compound.formula : `${coeff}${compound.formula}`;
        }).join(' + ').replace(/\+ (?=.*[A-Z].*→)/, ' → ');

        // Determine what changes were made
        const changes = [];
        result.compounds.forEach((compound, i) => {
            if (beforeCoeffs[i] !== afterCoeffs[i]) {
                changes.push({
                    compound: compound.formula,
                    before: beforeCoeffs[i],
                    after: afterCoeffs[i],
                    reason: `To balance ${element} atoms`
                });
            }
        });

        // Determine explanation
        let explanation = "";
        let title = "";
        let status = "";

        if (['Na', 'K', 'Ca', 'Mg', 'Fe', 'Cu', 'Zn', 'Al', 'Ag', 'Pb', 'Sn', 'Cr', 'Mn'].includes(element)) {
            title = `Balancing ${element} (Metal)`;
            status = "Step " + (index + 1) + " of " + balancingOrder.length;
            explanation = `Metals like ${element} usually appear in fewer compounds, making them easier to balance first. We need to adjust coefficients to make ${element} atoms equal on both sides.`;
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
            beforeEquation: beforeEquation,
            afterEquation: afterEquation,
            beforeLeftCompounds: beforeLeftCompounds,
            beforeRightCompounds: beforeRightCompounds,
            afterLeftCompounds: afterLeftCompounds,
            afterRightCompounds: afterRightCompounds,
            beforeLeftCount: beforeLeftCount,
            beforeRightCount: beforeRightCount,
            afterLeftCount: afterLeftCount,
            afterRightCount: afterRightCount,
            beforeBalanced: beforeLeftCount === beforeRightCount,
            afterBalanced: afterLeftCount === afterRightCount,
            changes: changes,
            explanation: explanation,
            // Fallback for old template
            currentEquation: beforeEquation,
            leftCompounds: beforeLeftCompounds,
            rightCompounds: beforeRightCompounds,
            leftCount: beforeLeftCount,
            rightCount: beforeRightCount,
            balanced: beforeLeftCount === beforeRightCount
        });

        // Update current coefficients for next iteration
        currentCoeffs = [...afterCoeffs];
    });

    return steps;
}