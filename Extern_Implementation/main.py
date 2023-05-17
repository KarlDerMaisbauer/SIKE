# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 14:11:42 2022

@author: David1
"""
from sidh_434 import*
import random as r
import math
#Parameters
#434
#xQ20 = 0x0000C7461738340EFCF09CE388F666EB38F7F3AFD42DC0B664D9F461F31AA2EDC6B4AB71BD42F4D7C058E13F64B237EF7DDD2ABC0DEB0C6C
#xQ21 = 0x000025DE37157F50D75D320DD0682AB4A67E471586FBC2D31AA32E6957FA2B2614C4CD40A1E27283EAAF4272AE517847197432E2D61C85F5
#yQ20 = 0x0001D407B70B01E4AEE172EDF491F4EF32144F03F5E054CEF9FDE5A35EFA3642A11817905ED0D4F193F31124264924A5F64EFE14B6EC97E5
#yQ21 = 0x0000E7DEC8C32F50A4E735A839DCDB89FE0763A184C525F7B7D0EBC0E84E9D83E9AC53A572A25D19E1464B509D97272AE761657B4765B3D6
#xP20 = 0x00003CCFC5E1F050030363E6920A0F7A4C6C71E63DE63A0E6475AF621995705F7C84500CB2BB61E950E19EAB8661D25C4A50ED279646CB48
#xP21 = 0x0001AD1C1CAE7840EDDA6D8A924520F60E573D3B9DFAC6D189941CB22326D284A8816CC4249410FE80D68047D823C97D705246F869E3EA50
#yP20 = 0x0001AB066B84949582E3F66688452B9255E72A017C45B148D719D9A63CDB7BE6F48C812E33B68161D5AB3A0A36906F04A6A6957E6F4FB2E0
#yP21 = 0x0000FD87F67EA576CE97FF65BF9F4F7688C4C752DCE9F8BD2B36AD66E04249AAF8337C01E6E4E1A844267BA1A1887B433729E1DD90C7DD2F
#xR20 = 0x0000F37AB34BA0CEAD94F43CDC50DE06AD19C67CE4928346E829CB92580DA84D7C36506A2516696BBE3AEB523AD7172A6D239513C5FD2516
#xR21 = 0x000196CA2ED06A657E90A73543F3902C208F410895B49CF84CD89BE9ED6E4EE7E8DF90B05F3FDB8BDFE489D1B3558E987013F9806036C5AC
#xQ30 = 0x00012E84D7652558E694BF84C1FBDAAF99B83B4266C32EC65B10457BCAF94C63EB063681E8B1E7398C0B241C19B9665FDB9E1406DA3D3846
#xQ31 = 0x00000000
#yQ30 = 0x00000000
#yQ31 = 0x0000EBAAA6C731271673BEECE467FD5ED9CC29AB564BDED7BDEAA86DD1E0FDDF399EDCC9B49C829EF53C7D7A35C3A0745D73C424FB4A5FD2
#xP30 = 0x00008664865EA7D816F03B31E223C26D406A2C6CD0C3D667466056AAE85895EC37368BFC009DFAFCB3D97E639F65E9E45F46573B0637B7A9
#xP31 = 0x00000000
#yP30 = 0x00006AE515593E73976091978DFBD70BDA0DD6BCAEEBFDD4FB1E748DDD9ED3FDCF679726C67A3B2CC12B39805B32B612E058A4280764443B
#yP31 = 0x00000000
#xR30 = 0x0001CD28597256D4FFE7E002E87870752A8F8A64A1CC78B5A2122074783F51B4FDE90E89C48ED91A8F4A0CCBACBFA7F51A89CE518A52B76C
#xR31 = 0x000147073290D78DD0CC8420B1188187D1A49DBFA24F26AAD46B2D9BB547DBB6F63A760ECB0C2B20BE52FB77BD2776C3D14BCBC404736AE4


#503
#xQ20 = 0x325CF6A8E2C6183A8B9932198039A7F965BA8587B67925D08D809DBF9A69DE1B621F7F134FA2DAB82FF5A2615F92CC71419FFFAAF86A290D604AB167616461
#xQ21 = 0x3E7B0494C8E60A8B72308AE09ED34845B34EA0911E356B77A11872CF7FEEFF745D98D0624097BC1AD7CD2ADF7FFC2C1AA5BA3C6684B964FA555A0715E57DB1
#yQ20 = 0x3A34654000BD4CB2612017BD5A1965A9F89FE11C55D517DF91B89B94F4F9C58B9A9DD056915573FEDC09CCD4997E82378759E00A2DE225CE04589D201FD754
#yQ21 = 0x19DEF0E8930E5123A22E346B1FFBD35EB01451647D8708A4835473B2539BD26806ED105A29F2D3F7EAA262426A965338782C5D20FBF478E4D1C8DBFC5B8294
#xP20 = 0x2ED31A03825FA14BC1D92C503C061D843223E611A92D7C5FBEC0F2C915EE7EEE73374DF6A1161EA00CDCB786155E21FD38220C3772CE670BC68274B851678
#xP21 = 0x1EE4E4E9448FBBAB4B5BAEF280A99B7BF86A1CE05D55BD603C3BA9D7C08FD8DE7968B49A78851FFBC6D0A17CB2FA1B57F3BABEF87720DD9A489B5581F915D2
#yP20 = 0x244D5F814B6253688138E317F24975E596B09BB15C6418E5295AAF73BA7F96EFED145DFAE1B21A8B7B121FEFA1B6E8B52F00478218589E604B97359B8A6E0F
#yP21 = 0x181CCC9F0CBE1390CC14149E8DE88EE79992DA32230DEDB25F04FADE07F242A9057366060CB59927DB6DC8B20E6B15747156E3C5300545E9674487AB393CA7
#xR20 = 0x3D24CF1F347F1DA54C1696442E6AFC192CEE5E320905E0EAB3C9D3FB595CA26C154F39427A0416A9F36337354CF1E6E5AEDD73DF80C710026D49550AC8CE9F
#xR21 = 0x6869EA28E4CEE05DCEE8B08ACD59775D03DAA0DC8B094C85156C212C23C72CB2AB2D2D90D46375AA6D66E58E44F8F219431D3006FDED7993F51649C029498
#xQ30 = 0x39014A74763076675D24CF3FA28318DAC75BCB04E54ADDC6494693F72EBB7DA7DC6A3BBCD188DAD5BECE9D6BB4ABDD05DB38C5FBE52D985DCAF74422C24D53
#xQ31 = 0x0
#yQ30 = 0x0
#yQ31 = 0x25512012C90A6869C4B29B9A757A03006BC7DF0BF7A2526A0713939FA48018AE3E249BD63699BEB3B8DEA215B7AE1B5A30FE371B64C5F1B0BF051A11D68E04
#xP30 = 0x32D03FD1E99ED0CB05C0707AF74617CBEA5AC6B75905B4B54B1B0C2D73697840155E7B1005EFB02B5D02797A8B66A5D258C76A3C9EF745CECE11E9A178BADF
#xP31 = 0x0
#yP30 = 0x2D810F828E3DC024D1BBBC7D6FA6E302CC5D458571763B7CCD0E4DBC9FA1163F0C1F8F4AE32A57F89DF8D2586D2A06E9FA30442B94A725266358C45236ADF3
#yP31 = 0x0
#xR30 = 0xC1465FD048FFB8BF2158ED57F0CFFF0C4D5A4397C7542D722567700FDBB8B2825CAB4B725764F5F528294B7F95C17D560E25660AD3D07AB011D95B2CB522
#xR31 = 0x288165466888BE1E78DB339034E2B8C7BDF0483BFA7AB943DFA05B2D1712317916690F5E713740E7C7D4838296E67357DC34E3460A95C330D5169721981758



#610
xQ20 = 0x25DA39EC90CDFB9BC0F772CDA52CB8B5A9F478D7AF8DBBA0AEB3E52432822DD88C38F4E3AEC0746E56149F1FE89707C77F8BA4134568629724F4A8E34B06BFE5C5E66E0867EC38B283798B8A
xQ21 = 0x2250E1959256AE502428338CB4715399551AEC78D8935B2DC73FCDCFBDB1A0118A2D3EF03489BA6F637B1C7FEE7E5F31340A1A537B76B5B736B4CDD284918918E8C986FC02741FB8C98F0A0ED
yQ20 = 0xA4FD5539025C0611E4B1CEC3C36F0D7590C035D3A25AD93022849CCEB3F67E4B1DBE988404290DD8B87B8D5E69ED3B0C5CDBCA248DC9D174CF762012CFE2D725CFD92057F2DBF8B04C7B12CC
yQ21 = 0x201C807BD738624E22B87554A2E053A46A9573BA863D4A9D309533E30B27BF7DD8137F5CE0F79C263D9D050541D69817A839085A76395F879315F6999E3441FC8FB3936DEE1BEF5B4E0E25096
xP20 = 0x1B368BC6019B46CD802129209B3E65B98BC64A92BC4DB2F9F3AC96B97A1B9C124DF549B528F18BEECB1666D27D47530435E84221272F3A97FB80527D8F8A359F8F1598D365744CA3070A5F26C
xP21 = 0x1459685DCA7112D1F6030DBC98F2C9CBB41617B6AD913E6523416CCBD8ED9C7841D97DF83092B9B3F2AF00D62E08DAD8FA743CBCCCC1782BE0186A3432D3C97C37CA16873BEDE01F0637C1AA2
yP20 = 0x1CD75CF512FFA9DF878EF495001A57ABC07FC7CE9BB488BB52DDCD7272D8A4FD17DD258ED3F844C862CF48803B9AC2668C7CB79C396128763B578080C30D14CA7EB709F98E3E682A391FB35A7
yP21 = 0x2001062A6289E4082CED884029207A1ACDEC525D7BC165A6CFF8BB469A8588950A416DBB924D2D673E3D6C32D232F6E6ADA62B37608F652C0B8628827B304BF1365D8211346207B24EFF09458
xR20 = 0x1B36A006D05F9E370D5078CCA54A16845B2BFF737C865368707C0DBBE9F5A62A9B9C79ADF11932A9FA4806210E25C92DB019CC146706DFBC7FA2638ECC4343C1E390426FAA7F2F07FDA163FB5
xR21 = 0x183C9ABF2297CA69699357F58FED92553436BBEBA2C3600D89522E7009D19EA5D6C18CFF993AA3AA33923ED93592B0637ED0B33ADF12388AE912BC4AE4749E2DF3C3292994DCF37747518A992
xQ30 = 0x14E647CB19B7EAAAC640A9C26B9C26DB7DEDA8FC9399F4F8CE620D2B2200480F4338755AE16D0E090F15EA1882166836A478C6E161C938E4EB8C2DD779B45FFDD17DCDF158AF48DE126B3A047
xQ31 = 0x0
yQ30 = 0x0
yQ31 = 0xE674067F5EA6DE85545C0A99E9E71E64FABFDC281D1E540FEDA47A56ED3ADCDDE1841083FABC7954B467C71AC6349B04974A7F9B688C5F735632FEB394146B0A080880069D8DA3324EDF795B
xP30 = 0x1587822E647707ED4313D3BE6A811A694FB201561111838A0816BFB5DEC625D23772DE48A26D78C04EEB26CA4A571C67CE4DC4C620282876B2F2FC2633CA548C3AB0C45CC991417A56F7FEFEB
xP31 = 0x0
yP30 = 0x14F295114B69D4769AC06DD07F051AD1114BCB7FA6B6EDE19F840969AFD56FD1F728907D3320A0309462A9444D24FE754666DB2470080951B31C2AC59704ABC7670C3C3A992C3C1629791F30
yP31 = 0x0
xR30 = 0x1DB73BC2DE666D24E59AF5E23B79251BA0D189629EF87E56C38778A448FACE312D08EDFB876C3FD45ECF3746D96E2CADBBA08B1A206C47DDD93137059E34C90E2E42E10F30F6E5F52DED74222
xR31 = 0x1B2C30180DAF5D91871555CE8EFEC76A4D521F877B754311228C7180A3E2318B4E7A00341FF99F34E35BF7A1053CA76FD77C0AFAE38E2091862AB4F1DD4C8D9C83DE37ACBA6646EDB4C238B48







Pa = point(fp2(xP20,xP21),fp2(1,0))
Qa = point(fp2(xQ20,xQ21),fp2(1,0))
Ra = point(fp2(xR20,xR21),fp2(1,0))
Pb = point(fp2(xP30,xP31),fp2(1,0))
Qb = point(fp2(xQ30,xQ31),fp2(1,0))
Rb = point(fp2(xR30,xR31),fp2(1,0))
#Initial Curve
a = fp2(6,0)
b = fp2(1,0)
AC = point(a,b)
###############################################################################
#Test : 3 ladder
#R = point()
#Enter value of sk :
#sk = 1
#R.ladder3pt(sk, Pa, Qa, Ra, AC)
#print("P + [sk]Q :",R)
"""
Can test other fonction too :
DBLADD, double, TPL
puissance2_k, TPLe
iso_2_curve, iso_3_curve, iso_4_curve
isogeny_2_point, isogeny_3_point, isogeny_4_point
"""



Pa = point(fp2(xP20,xP21),fp2(1,0))
Qa = point(fp2(xQ20,xQ21),fp2(1,0))
Ra = point(fp2(xR20,xR21),fp2(1,0))
Pb = point(fp2(xP30,xP31),fp2(1,0))
Qb = point(fp2(xQ30,xQ31),fp2(1,0))
Rb = point(fp2(xR30,xR31),fp2(1,0))


len = math.floor((math.log2(p) / 8))
#print(math.log2(p))
print(len)
shift = 0
#Alice_secret = 0x12452345235534652345678909876543234567899876543234567898765432345678998765432345678987654664745273466324008545bd06e2c2515f0074234654765865764c7545e14679e2a9e341b71efb2eb141f2507ed7ab3d1b58ba46e87ccd238e1f29625558ec2ae8944a19495cff74b0dc5166334873643c9869327b23c66b8b4566
Alice_secret = r.getrandbits(p.bit_length())
#Alice_secret = 0x1983f071f3b5e6abd907d045583239e5fc4d7769de42d1bcfd4a1cdf41c4d26937c27588d43ec7ccd238e1f29625558ec2ae8944a19495cff74b0dc5166334873643c9869648a3a1b4f5bb4e5
#while(Alice_secret.bit_length()//8 != 524):
#    Alice_secret = Alice_secret >> 1
Alice_secret = Alice_secret % (2**(e1))
#while Alice_secret.bit_length() != p.bit_length():
#    Alice_secret = r.getrandbits(p.bit_length())
#    Alice_secret = Alice_secret % p

print(Alice_secret.bit_length())
print(hex(p))
print(f'Secret alice \n{hex(Alice_secret)}')
x1 = fp2()
x2 = fp2()
x3 = fp2()

(x1, x2, x3) = isogen_2(Alice_secret, Pa.X, Qa.X, Ra.X, Pb.X, Qb.X, Rb.X)


print(" Alices publickey is")
print(f'{hex(x1.a0)} \n{hex(x1.a1)}')
print(f'{hex(x2.a0)} \n{hex(x2.a1)}')
print(f'{hex(x3.a0)} \n{hex(x3.a1)}')


x1_b = fp2()
x2_b = fp2()
x3_b = fp2()
#Bob_secret = 0x1436c6123245234543562bbd95a223456542345678909876543234567898765432345678987654323456789876543367547656757130a37c83e4583f2dba31431bd7b7519b500d25e43425234534255d324e6afb666b68079a41a7c4c91befd79f7fdcc2330ded7263109cf92e3352255a140e0f7666ef438d1190cde6
Bob_secret = r.getrandbits(p.bit_length())
#Bob_secret = 0x122d94421d217483a9d41294cde5350018b4d08179c0b2b2d9744d776c548b10ec963e40a6f06079a41a7c4c91befd79f7fdcc2330ded7263109cf92e3352255a140e0f7666ef438d393c4518
#while(Bob_secret.bit_length()//8 != 524):
#    Bob_secret = Bob_secret >> 1
s = len = math.floor((math.log2(3**e2)))
Bob_secret = Bob_secret % (2**(s))
#while Bob_secret.bit_length() != p.bit_length():
#    Bob_secret = r.getrandbits(p.bit_length())
#    Bob_secret = Bob_secret % p
print(f'Secret bob \n{hex(Bob_secret)}')
(x1_b, x2_b, x3_b) = isogen_3(Bob_secret, Pa.X, Qa.X, Ra.X, Pb.X, Qb.X, Rb.X)

print(f'bitlength p: {p.bit_length()}')
print(f'bitlength a: {Alice_secret.bit_length()}')
print(f'bitlength b: {Bob_secret.bit_length()}')


print("Bobs publickey is")
print(f'{hex(x1_b.a0)} \n{hex(x1_b.a1)}')
print(f'{hex(x2_b.a0)} \n{hex(x2_b.a1)}')
print(f'{hex(x3_b.a0)} \n{hex(x3_b.a1)}')



shared_alice = fp2()
shared_alice = isoex_2(x1_b, x2_b, x3_b, Alice_secret)
shared_alice.mod()
print("\n\n\nAlices' shared secret is:")
print(f'{hex(shared_alice.a0)} \n{hex(shared_alice.a1)}')


shared_bob = fp2()
shared_bob = isoex_3(x1, x2, x3, Bob_secret)
shared_bob.mod()
print("\nBobs' shared secret is:")
print(f'{hex(shared_bob.a0)} \n{hex(shared_bob.a1)}')


if shared_alice == shared_bob:
    print('####################')
    print('its the same picture')
    print('####################')
else:
    print('####################')
    print('not the same')
    print('####################')
