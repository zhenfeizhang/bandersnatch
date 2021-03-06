
# modulus
p=0x1ae3a4617c510eac63b05c06ca1493b1a22d9f300f5138f1ef3622fba094800170b5d44300000008508c00000000001
Fp=Zmod(p)

#####################################################
# Weierstrass curve: y² = x³ + A * x + B
#####################################################
# curve y^2 = x^3 + 1
WA=Fp(0)
WB=Fp(1)

#####################################################
# Montgomery curve: By² = x³ + A * x² + x
#####################################################
# root for x^3 + 1 = 0
alpha = -1
# s = 1 / (sqrt(3alpha^2 + a))
s = 1/(Fp(3).sqrt())

# MA = 3 * alpha * s
MA=Fp(228097355113300204138531148905234651262148041026195375645000724271212049151994375092458297304264351187709081232384)
# MB = s
MB=Fp(10189023633222963290707194929886294091415157242906428298294512798502806398782149227503530278436336312243746741931)

# #####################################################
# # Twised Edwards curve 1: a * x² + y² = 1 + d * x² * y²
# #####################################################
# a = (MA+2)/MB
TE1a=Fp(61134141799337779744243169579317764548490943457438569789767076791016838392692895365021181670618017873462480451583)
# b = (MA-2)/MB
TE1d=Fp(197530284213631314266409564115575768987902569297476090750117185875703629955647927409947706468955342250977841006588)

# #####################################################
# # Twised Edwards curve 2: a * x² + y² = 1 + d * x² * y²
# #####################################################
# a = -1
TE2a=Fp(-1)
# b = -TE1d/TE1a
TE2d=Fp(122268283598675559488486339158635529096981886914877139579534153582033676785385790730042363341236035746924960903179)


################################################################################
################################################################################
################################################################################
################################################################################

#####################################################
# Weierstrass curve generator
#####################################################
# obtained from arkworks code
Wx = Fp(81937999373150964239938255573465948239988671502647976594219695644855304257327692006745978603320413799295628339695)
Wy = Fp(241266749859715473739788878240585681733927191168601896383759122102112907357779751001206799952863815012735208165030)

assert(Wy^2 - Wx^3 - WA * Wx - WB == 0)

#####################################################
# Montgomery curve generator
#####################################################
# x = s * (x - alpha)
Mx = Fp(251803586774461569862800610331871502335378228972505599912537082323947581271784390797244487924068052270360793200630)
# y = s * y
My = Fp(77739247071951651095607889637653357561348174979132042929587539214321586851215673796661346812932566642719051699820)

assert(MB * My^2 == Mx^3+ MA * Mx^2 + Mx)

# #####################################################
# # Twised Edwards curve 1 generator
# #####################################################
# x = Mx/My
TE1x = Fp(82241236807150726090333472814441006963902378430536027612759193445733851062772474760677400112551677454953925168208)
# y = (Mx - 1)/(Mx+1)
TE1y = Fp(6177051365529633638563236407038680211609544222665285371549726196884440490905471891908272386851767077598415378235)

assert( TE1a * TE1x^2 + TE1y^2 == 1 + TE1d * TE1x^2 * TE1y^2 )


# #####################################################
# # Twised Edwards curve 2 generator
# #####################################################
beta = (-TE1a).sqrt()
# x = TE1x * sqrt(-TE1a)
TE2x = Fp(t)
# y = TE1y
TE2y = Fp(6177051365529633638563236407038680211609544222665285371549726196884440490905471891908272386851767077598415378235)

assert( TE2a * TE2x^2 + TE2y^2 == 1 + TE2d * TE2x^2 * TE2y^2 )
