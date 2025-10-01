# Chemical Equation Balancer

A web application that automatically balances chemical equations and provides step-by-step explanations of the balancing process.

## Features

- **Smart Input Parsing**: Accepts equations in various formats (H2 + O2 -> H2O, CH4 + O2 = CO2 + H2O)
- **Automatic Balancing**: Uses matrix-based linear algebra to find the smallest integer coefficients
- **Step-by-Step Explanations**: Shows detailed breakdown of the balancing process
- **Clean Interface**: Modern, responsive design with mobile support
- **Copy to Clipboard**: Easy sharing of balanced equations
- **Error Handling**: Clear feedback for invalid inputs

## How to Use

1. Open `index.html` in your web browser
2. Enter an unbalanced chemical equation in the input field
3. Click "Balance" or press Enter
4. View the balanced equation and step-by-step explanation

## Supported Input Formats

- Arrow notations: `->`, `=>`, `=`
- Spaces are optional: `H2+O2->H2O` or `H2 + O2 -> H2O`
- Existing coefficients are automatically removed and recalculated

## Examples

- Simple: `H2 + O2 -> H2O`
- Combustion: `CH4 + O2 -> CO2 + H2O`
- Complex: `C2H6 + O2 = CO2 + H2O`
- Oxidation: `Fe + O2 -> Fe2O3`

## Algorithm

The balancer uses a matrix-based approach:

1. **Parse** the equation into compounds and elements
2. **Build** a matrix representing the conservation of mass for each element
3. **Solve** the system of linear equations using Gaussian elimination
4. **Convert** to smallest positive integer coefficients
5. **Verify** the solution balances all elements

## Technical Details

- **Frontend**: HTML5, CSS (Tailwind), Vanilla JavaScript
- **Algorithm**: Linear algebra with Gaussian elimination
- **Performance**: Instant balancing (<1 second)
- **Compatibility**: Works in all modern browsers

## File Structure

```
├── index.html              # Main HTML file
├── equation-balancer.js     # Core balancing logic and UI
└── README.md               # This file
```

## Browser Support

- Chrome/Edge 60+
- Firefox 55+
- Safari 12+
- Mobile browsers

## License

MIT License - feel free to use and modify as needed.