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

    private readonly double _mass = 9.11e-31;
    private readonly double _hbar = 1.0545718e-34; // Reduced Planck's constant

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

    public Complex GaussianWavepacket(double x, double t, double x0, double sigma, double k)
    {
        // Gaussian wavepacket: ψ(x,t) = (2πσ²)^(-1/4) e^{-(x-x0-ħkt/m)²/(4σ²)} e^{i(kx - ωt)}
        var dispersion = sigma * Math.Sqrt(1 + (_hbar * t / (2 * _mass * sigma * sigma)) * (_hbar * t / (2 * _mass * sigma * sigma)));
        var norm = Math.Pow(2 * Math.PI * dispersion * dispersion, -0.25);
        
        var realExp = -Math.Pow(x - x0 - _hbar * _k * t / _mass, 2) / (4 * dispersion * dispersion);
        var phase = Complex.ImaginaryOne * (_k * x - _omega * t);
        
        return norm * Complex.Exp(realExp + phase);
    }

    public Complex Superposition(double x, double t, double k1, double k2)
    {
        // ψ = (1/√2) (e^{i(k1x - ω1t)} + e^{i(k2x - ω2t)})
        var omega1 = (_hbar * k1 * k1) / (2 * _mass);
        var omega2 = (_hbar * k2 * k2) / (2 * _mass);
        
        var psi1 = Complex.Exp(Complex.ImaginaryOne * (k1 * x - omega1 * t));
        var psi2 = Complex.Exp(Complex.ImaginaryOne * (k2 * x - omega2 * t));
        
        return (psi1 + psi2) * (1 / Math.Sqrt(2)); // Normalized
    }
}