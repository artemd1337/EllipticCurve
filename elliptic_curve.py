import gmpy2
import math


class EllipticCurve:
    def __init__(
            self,
            a: gmpy2.mpz,
            b: gmpy2.mpz,
            p: gmpy2.mpz
    ):
        assert a.bit_count() < 1024 and b.bit_count() < 1024 and p.bit_count() < 1024
        assert (gmpy2.mod(4 * gmpy2.powmod(a, 3, p), p) +
                gmpy2.mod(27 * gmpy2.powmod(b, 2, p), p) != 0)  # Сингулярная кривая 4a^3 + 27B^2 == 0
        assert p > 3
        self.a = gmpy2.mod(a, p)
        self.b = gmpy2.mod(b, p)
        self.p = p

    def __eq__(self, other):
        return self.a == other.a and self.b == other.b and self.p == other.p

    def calculate_order(self):
        if self.p <= 229:   # Если p достаточно мал, то return p + 1 + Sum(legendre(x^3+ax+b, p))
            return self.p + 1 + gmpy2.fsum(
                gmpy2.legendre(
                    gmpy2.mod(gmpy2.powmod(x, 3, self.p) + gmpy2.mul(self.a, x) + self.b, self.p),
                    self.p)
                for x in range(0, self.p)
            )

        state = gmpy2.random_state(hash(gmpy2.random_state()))
        g = gmpy2.mpz_random(state, self.p)  # генерируем случайный квадратичный невычет
        while gmpy2.legendre(g, self.p) != -1:
            g = gmpy2.mpz_random(state, self.p)
        W = math.ceil(math.pow(self.p, 1 / 4) * math.sqrt(2))  # Параметр giant step
        c = gmpy2.mod(gmpy2.mul(gmpy2.powmod(g, 2, self.p), self.a), self.p)  # Параметр a искаженной кривой
        d = gmpy2.mod(gmpy2.mul(gmpy2.powmod(g, 3, self.p), self.b), self.p)  # Параметр b искаженной кривой
        while True:
            x = gmpy2.mpz_random(state, self.p)  # Выбрать случайный x
            sigma = gmpy2.legendre(
                gmpy2.mod(gmpy2.powmod(x, 3, self.p) + gmpy2.mul(self.a, x) + self.b, self.p),
                self.p)
            #print(x, g, sigma)
            if sigma == 0:  # x^3+ax+b делится на p
                continue
            elif sigma == 1:  # x^3+ax+b - квадратичный вычет по модулю p -> Кривую не нужно искажать
                E = EllipticCurve(self.a, self.b, self.p)
            else:  # x^3+ax+b - квадратичный невычет -> исказить кривую и сделать допустимый x
                E = EllipticCurve(c, d, self.p)
                x = gmpy2.mod(gmpy2.mul(x, g), self.p)
            y = get_sqrt(gmpy2.mod(gmpy2.powmod(x, 3, E.p) + gmpy2.mul(E.a, x) + E.b, E.p), E.p)
            # print(f"x={x}, y={y}, E.a = {E.a}, E.b = {E.b}")
            P = EllipticCurvePoint(x=x, y=y, curve=E)

            # Рассчет пересечения двух списков
            A = [(P * (self.p + 1 + baby)) * x for baby in range(0, W)]
            """
            print("----")
            for baby in range(0, W):
                P_ = P.double_and_add(self.p + 1 + baby)
                Pn = P_.double_and_add(x)
                print(f"P={P}, P*({self.p}+1+{baby})={P_}, P*()*x={Pn}")
            print("----")
            """
            B = [(P * (gmpy2.mul(giant, W))) * x for giant in range(0, W + 1)]
            S = list(set(A) & set(B))
            """
            for elem in A:
                print(elem)
            print("---")
            for elem in B:
                print(elem)
            print("---")
            for elem in S:
                print(elem)
            print("---")
            """
            if len(S) != 1:
                continue
            s = S[0]
            baby_idx = A.index(s)
            giant_idx = B.index(s)
            print(f"baby_idx={baby_idx}, giant_idx={giant_idx}")
            t = baby_idx + gmpy2.mul(giant_idx, W)
            print(f"E.a = {E.a}, E.b = {E.b}")
            print(f"case1: t={t}, ord={self.p + 1 + gmpy2.mul(sigma, t)}, sigma={sigma}, (p + 1 + t)P = {P.ternary_mul(self.p + 1 + t)}")

            if not P.ternary_mul(self.p + 1 + t).is_inf:
                t = baby_idx - gmpy2.mul(giant_idx, W)
                print(f"case2: t={t}, ord={self.p + 1 + gmpy2.mul(sigma, t)}, sigma={sigma}, (p + 1 + t)P = {P.ternary_mul(self.p + 1 + t)}")
            return self.p + 1 + gmpy2.mul(sigma, t)

    def generate_point(self):
        state = gmpy2.random_state(hash(gmpy2.random_state()))
        while True:
            x = gmpy2.mpz_random(state, self.p)
            t = gmpy2.mod(gmpy2.mul(x, gmpy2.powmod(x, 2, self.p) + self.a) + self.b, self.p)
            if gmpy2.jacobi(t, self.p) == -1:
                continue
            else:
                return EllipticCurvePoint(
                    x=x,
                    y=get_sqrt(t, self.p),
                    curve=self,
                    is_inf=False
                )


def get_sqrt(x: gmpy2.mpz, p: gmpy2.mpz):
    """
    Нахождение квадратного корня по модулю числа p алгоритмом Тонелли-Шенкса
    Нотация из книги, страница 122
    """
    assert gmpy2.jacobi(x, p) == 1
    state = gmpy2.random_state(hash(gmpy2.random_state()))
    to_parse = p - 1
    S = 0
    for elem in gmpy2.digits(to_parse, 2)[::-1]:
        if elem == '1':
            break
        S += 1
    t = to_parse[S:]
    d = gmpy2.mpz_random(state, p)
    while gmpy2.jacobi(d, p) != -1:
        d = gmpy2.mpz_random(state, p)
    D = gmpy2.powmod(d, t, p)
    A = gmpy2.powmod(x, t, p)
    m = gmpy2.mpz(0)
    for i in range(1, S):
        if gmpy2.powmod(gmpy2.mul(A, gmpy2.powmod(D, m, p)), gmpy2.powmod(2, S - 1 - i, p), p) == gmpy2.mod(-1, p):
            m = gmpy2.add(m, gmpy2.powmod(2, i, p))
    result = gmpy2.mod(gmpy2.mul(
        gmpy2.powmod(x, gmpy2.divm(t + 1, 2, p), p),
        gmpy2.powmod(D, gmpy2.divm(m, 2, p), p)), p
    )
    return result


class EllipticCurvePoint:
    def __init__(
            self,
            x: gmpy2.mpz,
            y: gmpy2.mpz,
            curve: EllipticCurve,
            is_inf=False
    ):
        if is_inf:
            self.x = 0
            self.y = 0
        else:
            assert gmpy2.powmod(y, 2, curve.p) == gmpy2.mod(gmpy2.powmod(x, 3, curve.p) + curve.a * x + curve.b,
                                                            curve.p)
            self.x = x
            self.y = y

        self.curve = curve
        self.is_inf = is_inf

    def __add__(self, other: 'EllipticCurvePoint') -> 'EllipticCurvePoint':
        assert self.curve == other.curve
        if self.is_inf:
            return EllipticCurvePoint(
                x=other.x,
                y=other.y,
                curve=other.curve,
                is_inf=other.is_inf
            )

        if other.is_inf:
            return EllipticCurvePoint(
                x=self.x,
                y=self.y,
                curve=self.curve,
                is_inf=self.is_inf
            )

        if self == other:
            coeff = gmpy2.divm(3 * gmpy2.powmod(self.x, 2, self.curve.p) + self.curve.a, 2 * self.y, self.curve.p)
        else:
            if self.x == other.x:
                return EllipticCurvePoint(
                    x=gmpy2.mpz(0),
                    y=gmpy2.mpz(0),
                    curve=self.curve,
                    is_inf=True

                )
            else:
                coeff = gmpy2.divm(other.y - self.y, other.x - self.x, self.curve.p)

        result_x = gmpy2.mod(gmpy2.powmod(coeff, 2, self.curve.p) - self.x - other.x, self.curve.p)
        result_y = gmpy2.mod(coeff * (self.x - result_x) - self.y, self.curve.p)
        return EllipticCurvePoint(
            x=result_x,
            y=result_y,
            curve=self.curve
        )

    def __neg__(self) -> 'EllipticCurvePoint':
        return EllipticCurvePoint(
            x=self.x,
            y=gmpy2.mod(-self.y, self.curve.p),
            curve=self.curve
        )

    def __mul__(self, num: gmpy2.mpz) -> 'EllipticCurvePoint':
        return self.ternary_mul(num)

    def __eq__(self, other) -> bool:
        return self.x == other.x and self.y == other.y and self.curve == other.curve

    def __str__(self) -> str:
        return f"[{self.x}, {self.y}], if_inf: {self.is_inf}"

    def __hash__(self):
        return hash(self.x) + hash(self.y) + hash(self.is_inf) + hash(self.curve.a) + hash(self.curve.b) + hash(
            self.curve.p)

    def double_and_add(self, num: gmpy2.mpz) -> 'EllipticCurvePoint':
        if self.is_inf:
            return self
        if type(num) == int:
            num = gmpy2.mpz(num)
        result = EllipticCurvePoint(x=gmpy2.mpz(0), y=gmpy2.mpz(0), curve=self.curve, is_inf=True)
        for ind, bit in enumerate(num.digits(2)):
            result = result + result
            if bit == '1':
                result = result + self
        return result

    def ternary_mul(self, num: gmpy2.mpz):
        if self.is_inf:
            return self
        if type(num) == int:
            num = gmpy2.mpz(num)
        binary = num.digits(2)
        seq_len = 0
        ternary = []
        for bit in binary[::-1]:
            if bit == '0':
                if seq_len >= 2:
                    ternary.append(-1)
                    for i in range(seq_len - 1):
                        ternary.append(0)
                    seq_len = 1
                else:
                    for i in range(seq_len):
                        ternary.append(1)
                    ternary.append(0)
                    seq_len = 0
            else:
                seq_len += 1
        if seq_len >= 2:
            ternary.append(-1)
            for i in range(seq_len - 1):
                ternary.append(0)
            ternary.append(1)
        else:
            for i in range(seq_len):
                ternary.append(1)

        result = EllipticCurvePoint(x=0, y=0, curve=self.curve, is_inf=True)
        T = EllipticCurvePoint(x=self.x, y=self.y, curve=self.curve, is_inf=self.is_inf)
        for elem in ternary:
            if elem == 1:
                result = result + T
            elif elem == -1:
                result = result + (-T)
            T = T + T
        return result
