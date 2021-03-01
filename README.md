# FixedPointNumber
FixedPointNumber For Lua Runtime
64位定点数库，高位32位为整数部分, 最高位表示正负,低位32位为小数部分,整体用补码表示
表示范围-2^31 ~ 2^31-1,最小精度为1/(2^32) 
所有运算通过整型数的运算来实现,保证在所有机器上运行结果相同.
目前支持运算: 加减乘除,取最大最小,大小比较,开方,三角函数
目前支持的函数: __add,__sub,__mul,__div,__num,__tostring,__eq,__le, __lt,__le,equals,compare,min, max,clamp,tostring,tonumber,sqrt,ceil,floor,abs,sin, cos, tan,asin,acos, atan, atan2
常量 one,pi,zero, rad2deg, deg2rad, epsilon, piover2,pitimes2
