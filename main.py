import gmpy2
from elliptic_curve import EllipticCurve, EllipticCurvePoint
import unittest


class Test(unittest.TestCase):
    def test_sum(self) -> None:
        curve = EllipticCurve(a=gmpy2.mpz(0), b=gmpy2.mpz(2), p=gmpy2.mpz(7))
        P = EllipticCurvePoint(x=gmpy2.mpz(0), y=gmpy2.mpz(3), curve=curve)
        Q = EllipticCurvePoint(x=gmpy2.mpz(5), y=gmpy2.mpz(6), curve=curve)
        result = P + Q
        self.assertEqual(result.x, 6)
        self.assertEqual(result.y, 6)
        self.assertEqual(result.is_inf, False)

        P = EllipticCurvePoint(x=gmpy2.mpz(3), y=gmpy2.mpz(6), curve=curve)
        Q = EllipticCurvePoint(x=gmpy2.mpz(3), y=gmpy2.mpz(1), curve=curve)
        result = P + Q
        self.assertEqual(result.x, 0)
        self.assertEqual(result.y, 0)
        self.assertEqual(result.is_inf, True)

        curve = EllipticCurve(a=gmpy2.mpz(13), b=gmpy2.mpz(32), p=gmpy2.mpz(67))
        P = EllipticCurvePoint(x=gmpy2.mpz(7), y=gmpy2.mpz(8), curve=curve)
        Q = EllipticCurvePoint(x=gmpy2.mpz(64), y=gmpy2.mpz(10), curve=curve)
        result = P + Q
        self.assertEqual(result.x, 55)
        self.assertEqual(result.y, 15)
        self.assertEqual(result.is_inf, False)

        result = P + P
        self.assertEqual(result.x, 19)
        self.assertEqual(result.y, 6)
        self.assertEqual(result.is_inf, False)

        curve = EllipticCurve(
            a=gmpy2.mpz("123123679126794129398123981293123821939129321391293120123123123123"),
            b=gmpy2.mpz("1231231233284325982985987329749327432491283213812329839812932"),
            p=gmpy2.mpz("29933500047097369067937316004196040052873624302205671175238976482131293332251605769347243889197845099882939039929103427802794111"))
        P = EllipticCurvePoint(
            x=gmpy2.mpz("29244634081723867984438788561326585212900803192359039624106664878889385752078911721652892386416297514704264622266882233069590590"),
            y=gmpy2.mpz("3103463007522824575400885612643449683050049481367289110460009970054903906870734032856822506722779226673508674520703832292023951"),
            curve=curve,
            is_inf=False
        )
        Q = EllipticCurvePoint(
            x=gmpy2.mpz("9796373221671279900825788503809012723532997697622734107257883110269968394493381895266796306412784898753425436538244326923250928"),
            y=gmpy2.mpz("25801540760717384198196093129334839815036182246776913720794951456209662946864473005489364001000966660950994855007127484625120425"),
            curve=curve,
            is_inf=False
        )

        result = P + Q
        self.assertEqual(result.x, gmpy2.mpz("25184957971215549538784133806900255726296318962304994999627906554702859507316567119845830437405086520897168713605782708140686839"))
        self.assertEqual(result.y, gmpy2.mpz("11418188651670813883077556494707288656076989565581420853248143655109705181047795022477878321382564008989610294898802890161014306"))
        self.assertEqual(result.is_inf, False)

    def test_mul_double_and_add(self) -> None:
        curve = EllipticCurve(a=gmpy2.mpz(-2), b=gmpy2.mpz(7), p=gmpy2.mpz(19))
        P = EllipticCurvePoint(x=gmpy2.mpz(1), y=gmpy2.mpz(5), curve=curve)
        result = P.ternary_mul(gmpy2.mpz(15))
        self.assertEqual(result.x, 3)
        self.assertEqual(result.y, 3)
        self.assertEqual(result.is_inf, False)

        curve = EllipticCurve(a=gmpy2.mpz(13), b=gmpy2.mpz(32), p=gmpy2.mpz(67))
        P = EllipticCurvePoint(x=gmpy2.mpz(7), y=gmpy2.mpz(8), curve=curve)
        result = P.double_and_add(gmpy2.mpz(2))
        self.assertEqual(result.x, 19)
        self.assertEqual(result.y, 6)
        self.assertEqual(result.is_inf, False)

        curve = EllipticCurve(a=gmpy2.mpz(13), b=gmpy2.mpz(32), p=gmpy2.mpz(67))
        P = EllipticCurvePoint(x=gmpy2.mpz(7), y=gmpy2.mpz(8), curve=curve)
        result = P.double_and_add(gmpy2.mpz(13))
        self.assertEqual(result.x, 14)
        self.assertEqual(result.y, 12)
        self.assertEqual(result.is_inf, False)

        curve = EllipticCurve(a=gmpy2.mpz("123123679126794129398123981293123821939129321391293120123123123123"),
                              b=gmpy2.mpz("1231231233284325982985987329749327432491283213812329839812932"),
                              p=gmpy2.mpz(
                                  "29933500047097369067937316004196040052873624302205671175238976482131293332251605769347243889197845099882939039929103427802794111"))
        P = EllipticCurvePoint(
            x=gmpy2.mpz(
                "5583928004829285560570024109456618142543316055329415187993208279238079053938728416684206037947350498050330361104725123718913533"),
            y=gmpy2.mpz(
                "28670928789288006050981000648933177445923336416412088271696645015793621824621867492374076958598436491394481012062566436495787299"),
            curve=curve
        )
        result = P.double_and_add(13)
        self.assertEqual(
            result.x,
            gmpy2.mpz("24770582313186756921469767354027553514639009416746889224122575208947529756865280572245462183158554980035847621225819625268880844")
        )
        self.assertEqual(
            result.y,
            gmpy2.mpz("14407751359088024613596564170452036484226414110194546644742413667014372191725030214028430424181225174628186779026719618216247640")
        )
        self.assertEqual(result.is_inf, False)

        curve = EllipticCurve(
            a=gmpy2.mpz("123123679126794129398123981293123821939129321391293120123123123123"),
            b=gmpy2.mpz("1231231233284325982985987329749327432491283213812329839812932"),
            p=gmpy2.mpz("29933500047097369067937316004196040052873624302205671175238976482131293332251605769347243889197845099882939039929103427802794111"))
        P = EllipticCurvePoint(
            x=gmpy2.mpz("29244634081723867984438788561326585212900803192359039624106664878889385752078911721652892386416297514704264622266882233069590590"),
            y=gmpy2.mpz("3103463007522824575400885612643449683050049481367289110460009970054903906870734032856822506722779226673508674520703832292023951"),
            curve=curve,
            is_inf=False
        )
        result = P.double_and_add(gmpy2.mpz("9839842989872364387643768876436980430498498435688432709843798"))
        self.assertEqual(
            result.x,
            gmpy2.mpz("21235053704953676356096076306753008993950969725028202976671346849693517261122913181411988361840688070352589384751163857866960695")
        )
        self.assertEqual(
            result.y,
            gmpy2.mpz("28101692109526292993062488102117069016807791890719380932679029616771015757482614759261346977960341537861692898058486079562446642")
        )
        self.assertEqual(result.is_inf, False)
    def test_mul_ternary(self) -> None:
        curve = EllipticCurve(a=gmpy2.mpz(-2), b=gmpy2.mpz(7), p=gmpy2.mpz(19))
        P = EllipticCurvePoint(x=gmpy2.mpz(1), y=gmpy2.mpz(5), curve=curve)
        result = P.ternary_mul(gmpy2.mpz(15))
        self.assertEqual(result.x, 3)
        self.assertEqual(result.y, 3)
        self.assertEqual(result.is_inf, False)

        curve = EllipticCurve(a=gmpy2.mpz(13), b=gmpy2.mpz(32), p=gmpy2.mpz(67))
        P = EllipticCurvePoint(x=gmpy2.mpz(7), y=gmpy2.mpz(8), curve=curve)
        result = P.ternary_mul(gmpy2.mpz(2))
        self.assertEqual(result.x, 19)
        self.assertEqual(result.y, 6)
        self.assertEqual(result.is_inf, False)

        curve = EllipticCurve(a=gmpy2.mpz(13), b=gmpy2.mpz(32), p=gmpy2.mpz(67))
        P = EllipticCurvePoint(x=gmpy2.mpz(7), y=gmpy2.mpz(8), curve=curve)
        result = P.ternary_mul(gmpy2.mpz(13))
        self.assertEqual(result.x, 14)
        self.assertEqual(result.y, 12)
        self.assertEqual(result.is_inf, False)

        curve = EllipticCurve(a=gmpy2.mpz("123123679126794129398123981293123821939129321391293120123123123123"),
                              b=gmpy2.mpz("1231231233284325982985987329749327432491283213812329839812932"),
                              p=gmpy2.mpz(
                                  "29933500047097369067937316004196040052873624302205671175238976482131293332251605769347243889197845099882939039929103427802794111"))
        P = EllipticCurvePoint(
            x=gmpy2.mpz(
                "5583928004829285560570024109456618142543316055329415187993208279238079053938728416684206037947350498050330361104725123718913533"),
            y=gmpy2.mpz(
                "28670928789288006050981000648933177445923336416412088271696645015793621824621867492374076958598436491394481012062566436495787299"),
            curve=curve
        )
        result = P.ternary_mul(13)
        self.assertEqual(result.x, gmpy2.mpz(
            "24770582313186756921469767354027553514639009416746889224122575208947529756865280572245462183158554980035847621225819625268880844"))
        self.assertEqual(result.y, gmpy2.mpz(
            "14407751359088024613596564170452036484226414110194546644742413667014372191725030214028430424181225174628186779026719618216247640"))
        self.assertEqual(result.is_inf, False)

        curve = EllipticCurve(
            a=gmpy2.mpz("123123679126794129398123981293123821939129321391293120123123123123"),
            b=gmpy2.mpz("1231231233284325982985987329749327432491283213812329839812932"),
            p=gmpy2.mpz("29933500047097369067937316004196040052873624302205671175238976482131293332251605769347243889197845099882939039929103427802794111"))
        P = EllipticCurvePoint(
            x=gmpy2.mpz("29244634081723867984438788561326585212900803192359039624106664878889385752078911721652892386416297514704264622266882233069590590"),
            y=gmpy2.mpz("3103463007522824575400885612643449683050049481367289110460009970054903906870734032856822506722779226673508674520703832292023951"),
            curve=curve,
            is_inf=False
        )
        result = P.ternary_mul(gmpy2.mpz("9839842989872364387643768876436980430498498435688432709843798"))
        self.assertEqual(
            result.x,
            gmpy2.mpz("21235053704953676356096076306753008993950969725028202976671346849693517261122913181411988361840688070352589384751163857866960695")
        )
        self.assertEqual(
            result.y,
            gmpy2.mpz("28101692109526292993062488102117069016807791890719380932679029616771015757482614759261346977960341537861692898058486079562446642")
        )
        self.assertEqual(result.is_inf, False)

    def test_order(self) -> None:
        curve = EllipticCurve(a=gmpy2.mpz(-2), b=gmpy2.mpz(7), p=gmpy2.mpz(19))
        self.assertEqual(int(curve.calculate_order()), 21)

        curve = EllipticCurve(a=gmpy2.mpz(13), b=gmpy2.mpz(32), p=gmpy2.mpz(67))
        self.assertEqual(int(curve.calculate_order()), 56)

        curve = EllipticCurve(a=gmpy2.mpz(13), b=gmpy2.mpz(32), p=gmpy2.mpz(1009))
        self.assertEqual(int(curve.calculate_order()), 988)


def main():
    # Noting to do :(
    return 0


if __name__ == "__main__":
    main()
