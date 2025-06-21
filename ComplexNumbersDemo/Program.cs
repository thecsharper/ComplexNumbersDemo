using System.Numerics;

// Example parameters for an electron
var wavefunction = new QuantumWavefunction(k: 1.0, mass: 9.11e-31);

// Compute probability density at x=1.2 m, t=0.5 s 
Console.WriteLine($"Probability density: {wavefunction.ProbabilityDensity(1.2, 0.5)}"); 

public class QuantumWavefunction(double k, double mass, double hbar = 1.0545718e-34)
{
    // hbar Reduced Planck's constant (J·s)

    // Wavenumber (1/m)
    private readonly double _k = k;

    // Angular frequency (rad/s)
    private readonly double _omega = (hbar * k * k) / (2 * mass); 

    public Complex Compute(double x, double t)
    {
        if (t < 0) throw new ArgumentException("Time cannot be negative.");
   
        // ψ(x,t) = e^{i(kx - ωt)}, normalized for plane wave (non-physical normalization)
        return Complex.Exp(Complex.ImaginaryOne * (_k * x - _omega * t));
    }

    public double ProbabilityDensity(double x, double t)
    {
        var psi = Compute(x, t);
        
        return Complex.Abs(psi) * Complex.Abs(psi); // |ψ|²
    }
}