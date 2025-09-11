# Shamir's Secret Sharing - Polynomial Reconstruction

A Node.js implementation that reconstructs secrets from Shamir's Secret Sharing scheme using polynomial interpolation via the substitution method (Gaussian elimination).

## Overview

This tool solves Shamir's Secret Sharing problems by:
1. **Parsing encoded shares** from different number bases (binary, octal, decimal, hexadecimal, etc.)
2. **Reconstructing the polynomial** using the substitution method with exact BigInt arithmetic
3. **Extracting the secret** (constant term `a₀`) from the polynomial

The implementation uses exact fraction arithmetic to ensure precision and automatically excludes outlier points that don't fit the same polynomial.

## Features

- ✅ **Multi-base parsing** - Supports bases 2-36 with automatic conversion to decimal
- ✅ **Exact arithmetic** - Uses BigInt fractions for precise calculations
- ✅ **Outlier detection** - Automatically identifies and excludes inconsistent data points
- ✅ **Substitution method** - Implements polynomial interpolation via Vandermonde matrix solving
- ✅ **Clean output** - Returns only the secret constant in a simple format

## Mathematical Background

### Shamir's Secret Sharing
- A secret is encoded as the constant term (`a₀`) of a polynomial of degree `m = k-1`
- The polynomial is: `f(x) = aₘxᵐ + aₘ₋₁xᵐ⁻¹ + ... + a₁x + a₀`
- Shares are points `(x, y)` where `y = f(x)`
- Need at least `k` points to reconstruct the polynomial and recover `a₀`

### Substitution Method
For each point `(xᵢ, yᵢ)`, substitute `x = xᵢ` into the polynomial:
```
f(x₁) = aₘ(x₁)ᵐ + aₘ₋₁(x₁)ᵐ⁻¹ + ... + a₁(x₁) + a₀ = y₁
f(x₂) = aₘ(x₂)ᵐ + aₘ₋₁(x₂)ᵐ⁻¹ + ... + a₁(x₂) + a₀ = y₂
...
```

This creates a linear system of equations solved via Gaussian elimination to find the coefficients.

## Installation

```bash
# Clone or download the files
# Ensure you have Node.js installed (version 14+ recommended)
```

## Usage

### Basic Usage

```bash
node index.js sample.json
```

### Input Format

The JSON file should contain:
```json
{
  "keys": {
    "n": 10,  // Total number of shares
    "k": 7   // Minimum shares needed (threshold)
  },
  "1": {
    "base": "6",
    "value": "13444211440455345511"
  },
  "2": {
    "base": "15", 
    "value": "aed7015a346d635"
  }
  // ... more shares
}
```

### Output

The program outputs the reconstructed secret:
```
constant = 79836264049851
```

## Example

Given the sample data with k=7, n=10:
- **Input**: 10 shares in various bases (base-6, base-15, base-16, etc.)
- **Process**: Converts to decimal, excludes outliers (x=2, x=8), uses points {1,3,4,5,6,7,9}
- **Output**: `constant = 79836264049851`

## Algorithm Details

### 1. Base Conversion
Each share's value is converted from its specified base to decimal:
```javascript
// Example: "aed7015a346d635" in base-15 → decimal number
```

### 2. Outlier Exclusion
The algorithm automatically identifies consistent points that lie on the same polynomial, excluding outliers that would create inconsistent results.

### 3. Vandermonde Matrix Solution
Creates and solves the linear system:
```
[x₁ᵐ x₁ᵐ⁻¹ ... x₁ 1] [aₘ]   [y₁]
[x₂ᵐ x₂ᵐ⁻¹ ... x₂ 1] [aₘ₋₁] [y₂]
[  ⋮    ⋮   ⋱   ⋮ ⋮] [ ⋮ ] = [ ⋮]
[xₖᵐ xₖᵐ⁻¹ ... xₖ 1] [a₁]   [yₖ]
                       [a₀]
```

### 4. Secret Extraction
The secret is the constant term `a₀` from the reconstructed polynomial.

## Technical Implementation

- **Language**: Node.js (JavaScript)
- **Arithmetic**: BigInt for exact integer calculations
- **Fractions**: Custom Fraction class for precise rational arithmetic
- **Matrix Operations**: Gaussian elimination with partial pivoting
- **Base Parsing**: Supports bases 2-36 with digit validation

## File Structure

```
├── index.js      # Main implementation
├── sample.json   # Example input data
└── README.md     # This documentation
```

## Error Handling

The program handles various error conditions:
- Invalid JSON format
- Missing required fields (`keys.n`, `keys.k`)
- Invalid base digits
- Singular matrices (insufficient or inconsistent data)
- File I/O errors

## Performance Notes

- Uses exact arithmetic (no floating-point errors)
- Efficient for typical secret sharing parameters (k ≤ 20)
- Memory usage scales with polynomial degree and number precision
- Computation time is O(k³) for the matrix solution


