from typing import Callable, Final, Tuple
from typing_extensions import TypeAlias


IntervalEndPoints: TypeAlias = Tuple[float, float]
EulerNumberApproximatorFunction: TypeAlias = Callable[[None], float]


class TetrationIsNotDefinedWithArgOfNegativeTwoAndBelow(Exception):
    pass


precomputed_coeffs: Final[list[float]] = [
    1.0,
    1.091767351258322138,
    0.271483212901696469,
    0.212453248176258214,
    0.069540376139988952,
    0.044291952090474256,
    0.014736742096390039,
    0.008668781817225539,
    0.002796479398385586,
    0.001610631290584341,
    0.000489927231484419,
    0.000288181071154065,
    0.000080094612538551,
    0.000050291141793809,
    0.000012183790344901,
    0.000008665533667382,
    0.000001687782319318,
    0.000001493253248573,
    0.000000198760764204,
    0.000000260867356004,
    0.000000014709954143,
    0.000000046834497327,
    -0.000000001549241666,
    0.000000008741510781,
    -0.000000001125787310,
    0.000000001707959267,
]


def approx_e_helper_function() -> float:
    i = fact = result = 1
    _ = [result := result + 1 / (fact := fact * (i := i + 1)) for _ in 20 * "."]
    result += 1
    return result


def approx_base_e_tetration(
    arg: float,
    /,
    *,
    cutoff_point: float = 1.0,
    coeffs: list[float] = precomputed_coeffs,
    euler_number_approximator_function: EulerNumberApproximatorFunction = approx_e_helper_function,
) -> float:
    from functools import reduce

    if arg in [-1, 0]:
        return arg + 1
    if arg > cutoff_point:
        EULER_CONST: Final[float] = approx_e_helper_function()
        return EULER_CONST ** (
            approx_base_e_tetration(
                arg - 1,
                cutoff_point=cutoff_point,
                coeffs=coeffs,
                euler_number_approximator_function=euler_number_approximator_function,
            )
        )
    from math import log

    if -cutoff_point - 1 < arg < -cutoff_point:
        return log(
            approx_base_e_tetration(
                arg + 1,
                cutoff_point=cutoff_point,
                coeffs=coeffs,
                euler_number_approximator_function=euler_number_approximator_function,
            )
        )
    if -cutoff_point - 1 >= arg:
        raise TetrationIsNotDefinedWithArgOfNegativeTwoAndBelow
    solution: float = 0
    substituted_values: list[float] = [
        solution + coeff * (arg**index) for index, coeff in enumerate(coeffs)
    ]
    solution: float = reduce(lambda x, y: x + y, substituted_values)
    return solution


def plot_tetration_base_e_curve() -> None:
    import numpy as np
    import matplotlib.pyplot as plt

    STEP_LENGTH: Final[float] = 1e-4

    interval_end_points_of_plot: IntervalEndPoints = (-2 + 1e-5, 2.2)

    x = np.arange(
        interval_end_points_of_plot[0], interval_end_points_of_plot[1], STEP_LENGTH
    )
    y = list(map(lambda t: approx_base_e_tetration(t), x))

    _ = [plt.plot(x, y), plt.xlabel("x"), plt.ylabel("tet_e(x)"), plt.show()]


if __name__ == "__main__":
    """Approximate tetration with base e using MacLaurin series with precomputed coefficients and recursion"""
    plot_tetration_base_e_curve()
