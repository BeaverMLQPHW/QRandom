from qiskit import QuantumCircuit
from qiskit_aer import Aer
from qiskit_ibm_runtime import QiskitRuntimeService
from math import log as _log, exp as _exp, pi as _pi, e as _e, ceil as _ceil
from math import sqrt as _sqrt, acos as _acos, cos as _cos, sin as _sin
import requests


NV_MAGICCONST = 4 * _exp(-0.5)/_sqrt(2.0)
TWOPI = 2.0 * _pi
LOG4 = _log(4.0)
SG_MAGICCONST = 1.0 + _log(4.5)
BPF = 53        # Number of bits in a float
RECIP_BPF = 2**-BPF


class QRandom:
    def __init__(self, backend_name="qasm_simulator", use_ibm=None, token=None):
        self.backend_name = backend_name
        self.backend = Aer.get_backend(backend_name)
        self.gauss_next = None

        if use_ibm:
            QiskitRuntimeService.save_account(channel=backend_name, token=token, set_as_default=True, proxies=proxies)
            self.backend = QiskitRuntimeService()
        else:
            self.backend = Aer.get_backend(backend_name)

    def available_backends(self):
        '''Shows a list of available backends for the chosen provider.'''
        return self.backend.backends(proxies=proxies)

    def avaible_simulators(self):
        '''Shows a list of available Aer simulators.'''
        return Aer.backends()

    def _get_random_bits(self, n):
        """Generate n random bits using a quantum circuit."""
        qc = QuantumCircuit(n, n)
        qc.h(range(n))  # Apply Hadamard gate to all qubits
        qc.measure(range(n), range(n))  # Measure all qubits

        job = self.backend.run(qc, shots=1)
        result = job.result()
        counts = result.get_counts()
        bitstring = list(counts.keys())[0]

        return [int(bit) for bit in bitstring]

    def randint(self, a, b):
        """Generate a random integer between a and b (inclusive)."""
        difference = b - a
        n_bits = (difference).bit_length()
        while True:
            bits = self._get_random_bits(n_bits)
            rand_int = int("".join(map(str, bits)), 2)
            if 0 <= rand_int <= difference:
                return rand_int + a

    def random(self):
        """Generate a random float in the range [0.0, 1.0)."""
        n_bits = 53  # Number of bits for a double precision float
        bits = self._get_random_bits(n_bits)
        rand_int = int("".join(map(str, bits)), 2)
        return rand_int / (1 << n_bits)

    def choice(self, seq, weights=None):
        """Choose a random element from a non-empty sequence."""
        if not seq:
            raise IndexError("Cannot choose from an empty sequence.")

        if weights:
            random_weight = self.randint(0, sum(weights))
            index = weights.index(min(weights, key=lambda weight: abs(sum(weights[:weights.index(weight) + 1])- random_weight)))

        else:
            index = self.randint(0, len(seq) - 1)

        return seq[index]

    def shuffle(self, seq):
        """Shuffle the sequence in place."""
        for i in range(len(seq) - 1, 0, -1):
            j = self.randint(0, i)
            seq[i], seq[j] = seq[j], seq[i]

    def sample(self, population, k):
        """Choose k unique random elements from a population sequence."""
        len_population = len(population)

        if not 0 <= k <= len_population:
            raise ValueError("Sample larger than population or is negative")

        indices = []
        for _ in range(k):
            len_population -= 1
            random_index = self.randint(0, len_population)
            indices.append(population[random_index])
            population.pop(random_index)

        return indices

    def randrange(self, start, stop, step=1):
        """Choose a random item from range(start, stop[, step])."""
        difference = (stop - start - 1) // step

        n_bits = (difference).bit_length()

        while True:
            bits = self._get_random_bits(n_bits)
            rand_int = int("".join(map(str, bits)), 2)
            if 0 <= rand_int <= difference:
                return start + rand_int * step

    def uniform(self, a, b):
        '''Get a random number in the range [a, b) or [a, b] depending on rounding.'''
        return a + (b - a) * self.random()

    def triangular(self, low=0.0, high=1.0, mode=None):
        """Triangular distribution."""
        u = self.random()
        c = 0.5 if mode is None else (mode - low) / (high - low)

        if u > c:
            u = 1.0 - u
            c = 1.0 - c
            low, high = high, low

        return low + (high - low) * (u * c) ** 0.5

    def normalvariate(self, mu, sigma):
        """Normal distribution.

        mu is the mean, and sigma is the standard deviation.

        """

        random = self.random
        while 1:
            u1 = random()
            u2 = 1.0 - random()
            z = NV_MAGICCONST*(u1 - 0.5) / u2

            if (z ** 2 / 4.0) <= -_log(u2):
                break

        return mu + z * sigma

    def lognormalvariate(self, mu, sigma):
        """Log normal distribution.

        If you take the natural logarithm of this distribution, you'll get a
        normal distribution with mean mu and standard deviation sigma.
        mu can have any value, and sigma must be greater than zero.

        """
        return _exp(self.normalvariate(mu, sigma))

    def expovariate(self, lambd):
        """Exponential distribution.

        lambd is 1.0 divided by the desired mean.  It should be
        nonzero.  (The parameter would be called "lambda", but that is
        a reserved word in Python.)  Returned values range from 0 to
        positive infinity if lambd is positive, and from negative
        infinity to 0 if lambd is negative.

        """
        return -_log(1.0 - self.random()) / lambd

    def vonmisesvariate(self, mu, kappa):
        """Circular data distribution.

        mu is the mean angle, expressed in radians between 0 and 2*pi, and
        kappa is the concentration parameter, which must be greater than or
        equal to zero.  If kappa is equal to zero, this distribution reduces
        to a uniform random angle over the range 0 to 2*pi.

        """

        random = self.random
        if kappa <= 1e-6:
            return TWOPI * random()

        a = 1.0 + _sqrt(1.0 + 4.0 * kappa * kappa)
        b = (a - _sqrt(2.0 * a))/(2.0 * kappa)
        r = (1.0 + b * b)/(2.0 * b)

        while 1:
            u1 = random()

            z = _cos(_pi * u1)
            f = (1.0 + r * z)/(r + z)
            c = kappa * (r - f)

            u2 = random()

            if u2 < c * (2.0 - c) or u2 <= c * _exp(1.0 - c):
                break

        u3 = random()
        if u3 > 0.5:
            theta = (mu % TWOPI) + _acos(f)
        else:
            theta = (mu % TWOPI) - _acos(f)

        return theta

    def gammavariate(self, alpha, beta):
        """Gamma distribution.  Not the gamma function!

        Conditions on the parameters are alpha > 0 and beta > 0.

        The probability distribution function is:

                    x ** (alpha - 1) * math.exp(-x / beta)
          pdf(x) =  --------------------------------------
                      math.gamma(alpha) * beta ** alpha

        """

        if alpha <= 0.0 or beta <= 0.0:
            raise ValueError('gammavariate: alpha and beta must be > 0.0')

        random = self.random
        if alpha > 1.0:


            ainv = _sqrt(2.0 * alpha - 1.0)
            bbb = alpha - LOG4
            ccc = alpha + ainv

            while 1:
                u1 = random()
                if not 1e-7 < u1 < .9999999:
                    continue
                u2 = 1.0 - random()
                v = _log(u1/(1.0-u1))/ainv
                x = alpha*_exp(v)
                z = u1*u1*u2
                r = bbb+ccc*v-x
                if r + SG_MAGICCONST - 4.5*z >= 0.0 or r >= _log(z):
                    return x * beta

        elif alpha == 1.0:
            u = random()
            while u <= 1e-7:
                u = random()
            return -_log(u) * beta

        else:   # alpha is between 0 and 1 (exclusive)


            while 1:
                u = random()
                b = (_e + alpha)/_e
                p = b*u
                if p <= 1.0:
                    x = p ** (1.0/alpha)
                else:
                    x = -_log((b-p)/alpha)
                u1 = random()
                if p > 1.0:
                    if u1 <= x ** (alpha - 1.0):
                        break
                elif u1 <= _exp(-x):
                    break
            return x * beta

    def gauss(self, mu, sigma):
        """Gaussian distribution.

        mu is the mean, and sigma is the standard deviation.  This is
        slightly faster than the normalvariate() function.

        Not thread-safe without a lock around calls.

        """

        random = self.random
        z = self.gauss_next
        self.gauss_next = None
        if z is None:
            x2pi = random() * TWOPI
            g2rad = _sqrt(-2.0 * _log(1.0 - random()))
            z = _cos(x2pi) * g2rad
            self.gauss_next = _sin(x2pi) * g2rad

        return mu + z*sigma

    def betavariate(self, alpha, beta):
        """Beta distribution.

        Conditions on the parameters are alpha > 0 and beta > 0.
        Returned values range between 0 and 1.

        """

        y = self.gammavariate(alpha, 1.)
        if y == 0:
            return 0.0
        else:
            return y / (y + self.gammavariate(beta, 1.))

    def paretovariate(self, alpha):
        """Pareto distribution.  alpha is the shape parameter."""

        u = 1.0 - self.random()
        return 1.0 / pow(u, 1.0/alpha)

    def weibullvariate(self, alpha, beta):
        """Weibull distribution.

        alpha is the scale parameter and beta is the shape parameter.

        """

        u = 1.0 - self.random()
        return alpha * pow(-_log(u), 1.0/beta)

    def biased_bit(self, p):
        """Generate a 0 or 1 with probability p for 1 and (1-p) for 0 using a quantum circuit."""
        if not (0 <= p <= 1):
            raise ValueError("Probability p must be between 0 and 1.")

        theta = 2 * _acos(_sqrt(1 - p))

        qc = QuantumCircuit(1, 1)
        qc.ry(theta, 0)  # Apply Ry rotation to set bias
        qc.measure(0, 0)

        job = self.backend.run(qc, shots=1)
        result = job.result()
        counts = result.get_counts()
        return int(list(counts.keys())[0])

    def ANU_randint(self, start, end, key, length=1, dtype='uint8', blocksize=1):
        QRN_URL = "https://api.quantumnumbers.anu.edu.au/"

        params = {"length": length, "type": dtype, "size": blocksize}
        headers = {"x-api-key": key}
        num_type = int(''.join(i for i in dtype if i.isdigit()))

        response = requests.get(QRN_URL, headers=headers, params=params)

        if response.status_code == 200:
            js = response.json()
            dif = end - start
            if js["success"] == True:
                return [int(round(start + dif * i / num_type, 0)) for i in js["data"]]
            else:
                return js["message"]

        else:
            return response.text