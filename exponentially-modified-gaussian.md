# Exponentially modified Gaussian


## Literature

After reading through a bunch of papers that use Exponentially modified Gaussians,
I have not found one that is to my satisfaction.
Some highlights:

- [Grushka (1972). Characterization of exponentially modified Gaussian peaks in chromatography.] <br>
  **Problem:** Equations (1) and (2) cannot possibly be equivalent:
  - Equation (1) has dimensions of $A / \mathsf{time}$.
  - Equation (2) has dimensions of $A$.

- [Hanggi & Carr (1985). Errors in exponentially modified Gaussian equations in the literature.] <br>
  **Problem:**
  While its main result is correct (equations (1) through (4)),
  there are incorrect assertions in the steps taken to arrive at it
  (which is ironic, given that this is meant to be a paper that points out problems with other papers):
  - In (6), it specifies the unit area Gaussian function as being
    $1 / (\sigma \sqrt{2 \pi}) \exp (\dots)$ for $t > 0$, but zero for $t < 0$.
    This is incorrect; the Gaussian function should never be zero.
  - In (7), the expression for "exponential decay of unit area" has area $1/2$.

- [Kalambet et al. (2011). Reconstruction of chromatographic peaks using the exponentially modified Gaussian function.] <br>
  There is nothing wrong with this paper mathematically, but the choice of coefficient is strange.
  Rather than use a formula that has the area of the peak as a parameter,
  it gives an expression in terms of $h$, the height of the Gaussian before convolution with the exponential.
  This leads to a curve that has area $h \sigma \sqrt{2 \pi}$.

[Grushka (1972). Characterization of exponentially modified Gaussian peaks in chromatography.]:
  https://doi.org/10.1021/ac60319a011
[Hanggi & Carr (1985). Errors in exponentially modified Gaussian equations in the literature.]:
  https://doi.org/10.1021/ac00289a051
[Kalambet et al. (2011). Reconstruction of chromatographic peaks using the exponentially modified Gaussian function.]:
  https://doi.org/10.1002/cem.1343


## Derivation from first principles

The only way to be satisfied is to derive it yourself.

### Definition of exponentially modified Gaussian with unit area

Let $G(x)$ be the probability density for a normal distribution
with mean $\mu$
and standard deviation $\sigma$,

$$
  G(x) =
    \frac{1}{\sigma \sqrt{2 \pi}}
    \exp \left[
      -\frac{1}{2}
      \left(
        \frac{x - \mu}{\sigma}
      \right)^2
    \right].
$$

Let $E(x)$ be the probability density for an exponential distribution
with time scale $\tau$,

$$
  E(x) =
    \begin{cases}
      \frac{1}{\tau} \exp \left[ -\frac{x}{\tau} \right], & x \ge 0 \\
      0, & x < 0.
    \end{cases}
$$

Then the <i>exponentially modified Gaussian with unit area</i>
is defined as their convolution, i.e.

$$
\begin{align}
C(x)
  &= \int_{-\infty}^\infty G(x - x') E(x') \thinspace\mathrm{d}x' \\
  &=
    \int_0^\infty
      \frac{1}{\sigma \sqrt{2 \pi}}
      \exp \left[
        -\frac{1}{2}
        \left(
          \frac{x - x' - \mu}{\sigma}
        \right)^2
      \right]
      \frac{1}{\tau} \exp \left[ -\frac{x'}{\tau} \right]
    \thinspace\mathrm{d}x' \\
  &=
    \frac{1}{\tau \sigma \sqrt{2 \pi}}
    \int_0^\infty
      \exp \left[
        -\frac{1}{2}
        \left(
          \frac{x - \mu - x'}{\sigma}
        \right)^2
          -
        \frac{x'}{\tau}
      \right]
    \thinspace\mathrm{d}x'
\end{align}
$$

### Exp-erfc form

[Completing the square of the exponent], we get

$$
  -\frac{1}{2}
  \left(
    \frac{x - \mu - x'}{\sigma}
  \right)^2
    -
  \frac{x'}{\tau}
    =
  \frac{1}{2} \left( \frac{\sigma}{\tau} \right)^2
    - \frac{x - \mu}{\tau}
    - \frac{1}{2}
      \left(
        \frac{\sigma}{\tau} - \frac{x - \mu - x'}{\sigma}
      \right)^2.
$$

Therefore

$$
  C(x) =
    \frac{1}{\tau \sigma \sqrt{2 \pi}}
    \exp \left[
      \frac{1}{2} \left( \frac{\sigma}{\tau} \right)^2
      - \frac{x - \mu}{\tau}
    \right]
    \int_0^\infty
      \exp \left[
        -\frac{1}{2}
        \left(
          \frac{\sigma}{\tau} - \frac{x - \mu - x'}{\sigma}
        \right)^2
      \right]
    \thinspace\mathrm{d}x'.
$$

Next, we apply the change of variables
$z' = \frac{1}{\sqrt{2}} \left( \frac{\sigma}{\tau} - \frac{x - \mu - x'}{\sigma} \right)$:

- Exponent: $-\frac{1}{2} \left( \frac{\sigma}{\tau} - \frac{x - \mu - x'}{\sigma} \right)^2 = -z'^2$.
- Differential: $\mathrm{d}x' = \sigma\sqrt{2} \thinspace\mathrm{d}z'$.
- Lower bound: $x' = 0 \iff z' = \frac{1}{\sqrt{2}} \left( \frac{\sigma}{\tau} - \frac{x - \mu}{\sigma} \right)$.
- Upper bound: $x' = \infty \iff z' = \infty$.

Thus we have

$$
\DeclareMathOperator{\erfc}{erfc}
\begin{align}
  C(x)
  &=
    \frac{1}{\tau \sigma \sqrt{2 \pi}}
    \exp \left[
      \frac{1}{2} \left( \frac{\sigma}{\tau} \right)^2
      - \frac{x - \mu}{\tau}
    \right]
    \int_{z' = \frac{1}{\sqrt{2}} \left( \frac{\sigma}{\tau} - \frac{x - \mu}{\sigma} \right)}^\infty
      \exp \left[ -z'^2 \right]
    \sigma\sqrt{2} \thinspace\mathrm{d}z' \\
  &=
    \frac{1}{2 \tau}
    \exp \left[
      \frac{1}{2} \left( \frac{\sigma}{\tau} \right)^2
      - \frac{x - \mu}{\tau}
    \right]
    \frac{2}{\sqrt{\pi}}
    \int_{z' = \frac{1}{\sqrt{2}} \left( \frac{\sigma}{\tau} - \frac{x - \mu}{\sigma} \right)}^\infty
      \exp \left[ -z'^2 \right]
    \thinspace\mathrm{d}z' \\
  &=
    \frac{1}{2 \tau}
    \exp \left[
      \frac{1}{2} \left( \frac{\sigma}{\tau} \right)^2
      - \frac{x - \mu}{\tau}
    \right]
    \erfc \left[
      \frac{1}{\sqrt{2}} \left( \frac{\sigma}{\tau} - \frac{x - \mu}{\sigma} \right)
    \right].
\end{align}
$$

[Completing the square of the exponent]:
  https://tio.run/##y00syUjNTSzJTE78/z@1oiA/LzWvxCk1Lb8oVcFWQddQ30hBQ6NCQVchtxRIVAQUZeamairoKxRnpucmasYZwQX1SxJLrblgJjimlaQWAQ0A6werBclD1EON0wSJAFmoShDyCOv0YZZZ/wcK5JUoOCi4lebkBGfmFuRkplVGozvbVgHFGbHW/wE
  "Wolfram Language (Mathematica) – Try It Online"
