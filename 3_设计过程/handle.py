import os
from shapely.geometry import Polygon, LineString, Point, box
import math
import numpy as np
from stand import *
import plotly.graph_objs as go
import plotly
from matplotlib import pyplot as plt
import uuid

class Force():
	def __init__(self,num,degree,pos,moment=0):
		# num: 力的大小，degree: 力的方向（x轴夹角），pos: 力的位置Point类,moment: 弯矩
		self.num = num
		self.degree = degree
		self.pos = pos
		self.moment = moment
		self.y = self.num * math.sin(math.radians(self.degree))
		self.x = self.num * math.cos(math.radians(self.degree))
		self.vector = np.array([self.x,self.y]) # 力的向量

	def join(self,force):
		# 返回与其他力的结合力
		# force: 其他力，Force类
		# 返回值：结合力，Force类
		new_vector = self.vector + force.vector
		new_num = math.sqrt(np.sum(new_vector**2))
		moment = self.moment + force.moment #两个 力原有的弯矩之和
		moment += (force.pos.x - self.pos.x) * force.y # y方向分力导致的弯矩
		moment += (force.pos.y - self.pos.y) * force.x # x方向分力导致的弯矩
		degree = math.degrees(math.atan(new_vector[1]/new_vector[0]))
		return Force(new_num,degree,self.pos,moment=moment)
	def move_to(self,point):
		# 将力移动到某点 力的大小方向不变，弯矩、位置变化
		# point: 点，Point类
		# 返回值：力，Force类
		moment = self.moment + (self.pos.x - point.x) * self.y + (self.pos.y - point.y) * self.x
		return Force(self.num,self.degree,point,moment=moment)
	def print_data(self,desc = ''):
		if desc:
			print('\n',desc)
		print("""力的大小：{self.num}\n分力:[{self.x},{self.y}]\n力的方向(与x轴夹角)：{self.degree}\n力的位置：{self.pos.x},{self.pos.y}\n弯矩：{self.moment}""".format(self=self))
	def get_info(self):
		desc = 3
		return {
			'num':round(self.num,desc),
			"x":round(self.x,desc),
			"y":round(self.y,desc),
			# 'degree':self.degree,
			'pos':[round(self.pos.x,desc),round(self.pos.y,desc)]	
		}
# #Force 类测试
# a = Force(2,45,Point(0,0),3)
# b = Force(3,-90,Point(1,0))
# c = a.join(b)
# c.print_data()
# d = Force(1,-90,Point(1,0))
# d.move_to(Point(2,1)).print_data()
class Dangtuqiang():
	def __init__(self, data):
		# read data
		self.design = data["design"]
		self.tu = data["tu"]
		self.parms = data["parms"]
		self.phi = data["parms"]["phi"]
		self.accuracy = data["parms"]["accuracy"]
		self.fake_wall_line = LineString(
			[data["design"]["coord_qz"], data["design"]["coord_qd"]])  # 假想墙背线
		self.po_xian = LineString(data["tu"]["po"])
		self.qz = data["design"]["coord_qz"]
		self.qd = data["design"]["coord_qd"]
		self.gamma_wall = data["parms"]["gamma_wall"]
		self.gamma_soil = data["parms"]["gamma_soil"]
		#.eta = math.degrees(math.atan(self.parms["C_z"]*self.parms["k_H"])) # 地震角
		self.C_z = data["parms"]["C_z"]
		self.k_H = data["parms"]["k_H"]
		self.base_type = data["parms"]["base_type"]
		self.gamma_0 = data["parms"]["gamma_0"]
		self.i = 0
		# init
		self.alpha, self.theta = 0, 0
		self.dec = 3  # 保留小数点后几位
		self.rho = round(math.degrees(
			math.atan(data["design"]["B3"]/data["design"]["H"])),self.dec)
		
	def sin(self, degree):
		# 接受degree，转化为弧度计算sin
		return math.sin(math.radians(degree))

	def cos(self, degree):
		return math.cos(math.radians(degree))

	def tan(self, degree):
		return math.tan(math.radians(degree))

	def show(self, coords):
		# 根据 coords 显示图形
		x, y = [], []
		for i in coords:
			x.append(i[0])
			y.append(i[1])
		plt.plot(x, y)
	def show_po(self):
		# 调试 显示 原点到交线楔形体的面积
		self.show(self.tu["po"])
		for i in self.tu["tuzhu"]:
			self.show(i["coords"])
		# print("交点",jd1,jd2,list(jx.coords))
		# print("墙踵",qz)
		#self.show([self.qz]+list(self.jx.coords)+[self.qz])
		self.show(self.design["coords"])
	def show_shape(self,shape):
		# 传入shapely对象
		self.show(shape.exterior.coords)
	def cal_stress(self,coords):
		# 计算应力,输入土体高度，输出应力，《铁道工程》p291
		return [[self.lambda_a * self.gamma_soil * i[0],i[1]] for i in coords]
	
	def to_round(self,list_or_num):
		# 将数字或者数组转化为保留小数点后几位的数字
		if isinstance(list_or_num,list):
			return [round(i,self.dec) for i in list_or_num]
		else:
			return round(list_or_num,self.dec)
	def cal_jx(self,qz,alpha,theta):
		#返回破裂面之间的交线
		x1 = qz[0] - (100-qz[1]) * self.tan(alpha) #无限远处（100远）的第二破裂面x坐标
		x2 = qz[0] + (100-qz[1]) * self.tan(theta) #无限远处（100m远）的第1破裂面x坐标
		jx = Polygon([qz, [x1, 100], [x2, 100]]).intersection(self.po_xian)  # 中间交线，可能为空
		return jx
	def get_centroid(self,polygons):
		# 获取多个多边形的质心
		# 输入多边形列表
		x = 0
		y = 0
		total_area = 0
		for i in polygons:
			area = i.area
			total_area += area
			x += i.centroid.x*area
			y += i.centroid.y*area
		return Point(x/total_area,y/total_area)

	def cal_area(self, alpha, theta):
		# 给定alpha,theta，计算楔形体的面积

		# 过qz和qd点的直线与坡线的交点
		qz = self.qz #墙踵坐标，作为原点
		jx = self.cal_jx(qz,alpha,theta) #交线
		# # 调试 显示 原点到交线楔形体的面积
		# self.show(self.tu["po"])
		# for i in self.tu["tuzhu"]:
		# 	self.show(i["coords"])
		# # print("交点",jd1,jd2,list(jx.coords))
		# # print("墙踵",qz)
		# self.show([qz]+list(jx.coords)+[qz])
		# self.show_po()
		# plt.savefig(f"./{self.i}.png")
		# self.i += 1
		# #plt.show()
		# plt.close()


		#计算填土部分的多边形面积 qz->交线->qz
		xxt = Polygon([qz]+list(jx.coords)+[qz])

		# 计算土柱中属于第一二破裂面的部分的面积及属于稳定自重的部分
		before_tuzhu = []  # 土柱中属于自重的部分，保存成一个多边形，方便计算重心
		for i in self.tu["tuzhu"]:
			# 遍历每个土柱，判断土柱是否在交线上
			dx_tuzhu = LineString([i["start"],i["end"]]) #土柱底部线段	
			jx_tuzhu = dx_tuzhu.intersection(jx) #土柱底部线段与交线相交部分
			if jx_tuzhu.is_empty:
				#未相交，可能在左右两侧
				if i["end"][0]<= list(jx.coords)[0][0]:
					#土柱在左侧
					before_tuzhu += [Polygon(i["coords"])]
				else:
					#右侧不管
					pass
			else:
				#相交
				jx_tuzhu_start =list(jx_tuzhu.coords)[0]
				jx_tuzhu_end =list(jx_tuzhu.coords)[-1]
				xxt = xxt.union(Polygon([jx_tuzhu_start,jx_tuzhu_end,[jx_tuzhu_end[0],jx_tuzhu_end[1]+i["height"]],[jx_tuzhu_start[0],jx_tuzhu_start[1]+i["height"]],jx_tuzhu_start]))
				if i["start"][0] < jx_tuzhu_start[0]:
					#土柱有一部分在左侧
					before_tuzhu += [Polygon([i["start"],jx_tuzhu_start,[jx_tuzhu_start[0],jx_tuzhu_start[1]+i["height"]],[i["start"][0],i["start"][1]+i["height"]],i["start"]])]
					

		#before_tuzhu_area,#土柱中属于自重的面积
		# 楔形体全部面积为楔形体面积加土柱中属于楔形体的面积
		result = (before_tuzhu,xxt)
		return result
	
	def cal_force_pos(self,alpha,theta):
		# 输入alpha，theta，汇出土压力图形，计算作用点位置
		coords = [[0,0]]
		self.lambda_a = (self.tan(self.theta)-self.tan(self.alpha))*self.cos(self.theta+self.phi)/self.sin(self.theta+2*self.phi+self.alpha)
		# 计算楔形体中的土压力应力
		self.jx = self.cal_jx(self.qz,alpha,theta)
		jx_start = list(list(self.jx.coords)[0])
		jx_end = list(list(self.jx.coords)[-1])
		self.jx_start = jx_start
		self.jx_end = jx_end
		self.second_plm_line = LineString([self.qz,jx_start])
		self.first_plm_line = LineString([self.qz,jx_end])
		x1 =jx_start[0] + (100-jx_start[1])*self.tan(theta)
		pol = Polygon([self.qz,jx_start,[x1,100],jx_end,self.qz])
		tmp_jx = pol.intersection(self.po_xian)
		coords.append([0,jx_start[1]])
		for point in list(tmp_jx.coords)[1:-1]:
			#过point做夹角为(90-theta)的直线，交fake_wall_line于jd
			x1 = point[0] - point[1]*self.tan(theta)
			jd = LineString([[x1,0],point]).intersection(self.second_plm_line)
			# print("绘制交点",point,jd)
			# self.show([[x1,0],point])
			coords.append([point[1]-jd.y,jd.y])
		coords += [[jx_end[1],0],[0,0]]
		coords = self.cal_stress(coords)
		pol1 = Polygon(coords) #由填土构成的土压力应力图
		for i in self.tu["tuzhu"]:
			start,end,height = i["start"],i["end"],i["height"]
			if start[0]<jx_start[0]<end[0]:
				#土柱在左侧中间
				start = jx_start
			elif start[0]<jx_end[0]<end[0]:
				#土柱在右侧中间
				end = jx_end
			elif jx_start[0] <= start[0] and end[0] <= jx_end[0]:
				#土柱在楔形体内部
				pass
			else:
				#土柱在楔形体外部
				continue
			#计算土柱中构成的土压力应力图
			jd1 = LineString([start,[start[0]-start[1]*self.tan(theta),0]]).intersection(self.second_plm_line)
			jd2 = LineString([end,[end[0]-end[1]*self.tan(theta),0]]).intersection(self.second_plm_line)
			coords2 = [[0,jd1.y],[start[1]-jd1.y+height,jd1.y],[end[1]-jd2.y+height,jd2.y],[0,jd2.y]]
			coords2 = self.cal_stress(coords2)
			pol1 = pol1.union(Polygon(coords2))
		# 混凝土夹角
		# H3 = self.design["H3"]
		# jd3 = LineString([[self.qz[0]-H3*self.tan(theta),0],[self.qz[0],H3]]).intersection(self.second_plm_line)
		# pol3 = Polygon([[0,0],[0,H3],[(H3-jd3.y)*self.gamma_wall/self.gamma_soil,jd3.y],[0,0]])
		#centroid = self.get_centroid([pol1,pol3]) #3.1942848678412874 差别不大
		centroid = pol1.centroid
 		
		# # 应力图展示
		plt.axis("off")
		plt.rcParams['font.family'] = 'SimHei'
		plt.cla()
		plt.subplot(221)
		plt.title('挡土墙截面图')
		ax1 = plt.gca()
		ax1.set_ylim([0,20])
		ax1.set_xlim([0,20])
		self.show_po()
		plt.subplot(222)
		plt.title('土压力分布及重心位置')
		ax1 = plt.gca()
		ax1.set_ylim([0,20])
		ax1.set_xlim([0,50])
		self.show(list(pol1.exterior.coords))
		plt.scatter(centroid.x,centroid.y)
		plt.text(centroid.x+1,centroid.y-1,'重心')
		plt.margins(0, 0)
		# plt.show()
		self.stress_path =f"/tmp/{uuid.uuid4().hex}.png"
		plt.savefig("./html"+self.stress_path,dpi=200,bbox_inches='tight')
		plt.close()

		#计算土压力作用点`	`
		return Point(self.qz[0]-self.tan(alpha)*centroid.y,centroid.y)

	def cal_E_a(self):
		# 土压力、楔体重力
		tmps = []
		Z = []
		theta_range = np.arange(0, 90-self.phi, self.accuracy)
		if theta_range[-1] != 90-self.phi:
			theta_range = np.append(theta_range, 90-self.phi)
		alpha_range = np.arange(self.accuracy, self.rho, self.accuracy) #跳过0
		if alpha_range[-1] != self.rho:
			alpha_range = np.append(alpha_range, self.rho)
		theta_range = list(theta_range)
		alpha_range = list(alpha_range)
		##############################
		# alpha_range = [15.331]
		# theta_range = [25.812]
		##############################
		for theta in theta_range:
			z_s = []
			for alpha in alpha_range:
				alpha = round(alpha, self.dec)
				# print(theta, alpha)
				# 计算楔形体的重力
				self.before_tuzhu,xxt = self.cal_area(alpha, theta) #土柱中属于自重的面积，楔形体的面积
				W1 = Force(xxt.area* self.gamma_soil,-90,xxt.centroid)
				# 计算土压力
				E = self.sin(90-theta-self.phi)/self.sin(alpha+theta+2*self.phi)*W1.num
				Ex = E*self.sin(90-alpha-self.phi)
				Ey =  E*self.cos(90-alpha-self.phi)
				tmps.append((theta,alpha,Ex,Ey,E,W1,xxt.area))
				z_s.append(Ex)
			Z.append(z_s)

		tmp_max = max(tmps, key=lambda x: x[2])
		self.theta,self.alpha,self.E_x,self.E_y,self.E,self.W1,self.area_W1 = tmp_max
		self.alpha,self.theta = float(self.alpha),float(self.theta)
		self.have_2_plm = 0 if self.alpha == self.rho else 1 #是否有第二破裂面
		# 此时已经求出了土压力大小，方向已知
		# print('max', tmp_max)
		# # 2d 图
		# c = [[x[1], x[0]*self.sin(90-x[2]-self.phi)] for x in tmps]
		# self.show(c)
		# plt.show()

		#3d图
		# fig = plt.figure()  # 定义新的三维坐标轴
		# # ax3 = plt.axes(projection='3d')
		# # # 作图
		# # ax3.plot_surface(X, Y, Z, cmap='rainbow')
		# # self.png_E_3d = self.get_png(plt)
		# # # plt.show()
		# # plt.close()
		# # 定义三维数据

		X, Y = np.meshgrid(alpha_range, theta_range)
		Z = np.array(Z)
		data=[go.Surface(x=X, y=Y, z=Z)]
		fig = go.Figure(data)
		fig.update_layout(
			title_text='土压力随破裂角变化图',
			autosize=True,	
			scene=dict(
				xaxis=dict(title = "第二破裂角"),
				yaxis=dict(title = '第一破裂角'),
				zaxis=dict(title = '主动土压力水平分力'),
				annotations=[
				dict(
					x= self.alpha,
					y= self.theta,
					z=self.E_x,
					text="主动土压力Ea",
					textangle=0,
					ax=0,
					ay=-75,
					font=dict(
						color="black",
						size=12
					)
				)]
			),
		)
		self.png_E_3d = plotly.io.to_json(fig)
		# 求作用点,画应力图
		pos = self.cal_force_pos(self.alpha,self.theta)
		
		# 土压力解算完毕，大小为self.E,方向为 self.alpha+self.fhi,作用点为pos
		return Force(self.E,self.alpha+self.phi,pos)

	def cal_W(self):
		# 计算挡墙自重，包括混凝土墙身、第二破裂面与真实墙背间的土重、属于自重的土柱
		# 1.混凝土墙身,T形截面重力W1_t,扶壁重力算入土重
		pol1 = Polygon(self.design["coords"])
		# 2.第二破裂面与真实墙背间的土重
		H3 =  self.design["H3"]
		x1 = self.qz[0]-self.tan(self.alpha)*H3 #第二破裂面与墙踵上边缘的交点
		left_po_xian = box(0,0,self.jx_start[0],self.jx_start[1]+100).intersection(self.po_xian)
		pol2 = Polygon([[x1,H3],[self.qz[0]-self.design["B3"],H3]]+list(left_po_xian.coords)+[[x1,H3]])
		# 3.属于自重的土柱
		pol3 = self.before_tuzhu
		
		#绘图 重力部分图形
		# self.show_shape(W2)
		# self.show_shape(W1)
		# for i in W3:
		# 	self.show_shape(i)
		# plt.show()

		#总的重力 W,与理正误差为千分之1，忽略不计
		for i in pol3:
			pol2 = pol2.union(i)
		self.G = Force(pol1.area*self.gamma_wall,-90,pol1.centroid)
		self.W = Force(pol2.area*self.gamma_soil,-90,pol2.centroid)
		self.area_G = pol1.area
		self.area_W = pol2.area
		return self.G,self.W #G为挡墙自重，W实际墙背与第二破裂面之间的土重 W1为破裂体的重力

	def cal_exforce(self):
		# 计算外力
		# 1.土压力
		self.E_a = self.cal_E_a()
		# self.E_a.print_data("E_a")
		# 2.自重
		self.G,self.W = self.cal_W() #G为挡墙自重，W实际墙背与第二破裂面之间的土重
		#墙趾板土重力
		W_toe = Polygon([[0,self.design["H2"]],[self.design["B2"],self.design["Hzd"]],[self.design["B2"],self.design["cover_h"]],[0,self.design["cover_h"]],[0,self.design["H2"]]])
		self.W_toe = Force(W_toe.area*self.gamma_soil,-90,W_toe.centroid)

		# 3.地震时土压力
		# 第一二破裂面之间的楔形体的地震力
		self.W1_p = Force(self.W1.num * self.C_z *self.k_H,180,self.W1.pos)
		# 第二破裂面与真实墙背间的地震力
		self.W_p = Force(self.W.num * self.C_z *self.k_H,180,self.W.pos)
		# 挡墙自身受到的地震力
		self.G_p = Force(self.G.num * self.C_z *self.k_H,180,self.G.pos)

	def check_k_c(self,data,gamma_E1,gamma_G,gamma_E2):
		#滑动稳定性检测 《极限状态法》p91
		W,E_y,E_x,F_hE,f_p,Z_x,Z_y,Z_hE,Z_W = data
		S_ddst = gamma_E1 * E_x + F_hE #不平衡作用效应设计值
		S_dstb = (gamma_G * W +gamma_E2 *E_y)*f_p #平衡作用设计值
		return self.to_round([S_ddst*self.gamma_0,S_dstb])
	def check_k_0(self,data,gamma_E1,gamma_G,gamma_E2):
		#倾覆稳定性检测 《极限状态法》p93
		W,E_y,E_x,F_hE,f_p,Z_x,Z_y,Z_hE,Z_W = data
		S_ddst = gamma_E1 * E_x * Z_x  + F_hE * Z_hE  #不平衡作用效应设计值
		S_dstb = gamma_G * W * Z_W + gamma_E2 * E_y * Z_y #平衡作用设计值
		return self.to_round([S_ddst*self.gamma_0,S_dstb])
		
	def check_e(self,data,limit):
		W,E_y,E_x,F_hE,f_p,Z_x,Z_y,Z_hE,Z_W = data
		self.B = self.qz[0] #基底宽度
		#稳定力系
		self.sum_M_y =  self.W_toe.num * self.W_toe.pos.x + \
						self.G.num * self.G.pos.x + \
						self.W.num * self.W.pos.x + \
						abs(self.E_a.y) * self.E_a.pos.x
		#倾覆力系力矩
		self.sum_M_0 =  abs(self.E_a.x) * self.E_a.pos.y + \
						self.G_p.num * self.G_p.pos.y + \
						self.W_p.num * self.W_p.pos.y + \
						self.W1_p.num * self.E_a.pos.y #此处为简化

		self.sum_N =  W + abs(self.E_a.y) + self.W.num +self.W_toe.num #作用在基底上的竖直力 ：挡墙自重，主动土压力竖向分力，实际墙背与第二破裂面之间的土重
		self.c = (self.sum_M_y - self.sum_M_0)/self.sum_N
		self.e = self.B/2 - self.c
		limit =  self.B * limit #偏心距限定值 ，《极限状态法》p95
		return self.to_round([self.e,limit])

	def check_sigma(self,sigma,gamma_toe,gamma_heel,gamma_a):
		#基底应力检测
		# sigma 地基承载力特征值
		if abs(self.e) <= self.B/6:
			sigma_1k = self.sum_N/self.B*(1+6*self.e/self.B)
			sigma_2k = self.sum_N/self.B*(1-6*self.e/self.B)
		elif self.e > self.B/6:
			sigma_1k = 2*self.sum_N/3/self.c
			sigma_2k = 0
		else:
			sigma_1k = 0
			sigma_2k = 2*self.sum_N/3/(self.B - self.c)
		sigma_pk = (sigma_1k + sigma_2k)/2
		return {
				"sigma_1k":self.to_round([sigma_1k,gamma_toe * sigma]),
				"sigma_2k":self.to_round([sigma_2k,gamma_heel * sigma]),
				"sigma_pk":self.to_round([sigma_pk,gamma_a * sigma])
			}

	def check(self,type = 'I'):
		# 检测稳定性
		# 输入type为组合类型，I为永久+主可变，IV 为永久+主可变+地震
		W = self.G.num
		E_y = self.E_a.y + self.W.num + self.W_toe.num
		E_x = self.E_a.x
		F_hE = self.G_p.num #无地震时G是0
		if type == "IV":
			#验算地震时
			E_x = self.E_a.x + self.W_p.num + self.W1_p.num
		f_p = 1.5 *self.parms["f"]
		Z_x = self.E_a.pos.y
		Z_y = self.E_a.pos.x
		Z_hE = self.G_p.pos.y
		Z_W = self.G.pos.x
		data = (W,E_y,E_x,F_hE,f_p,Z_x,Z_y,Z_hE,Z_W)
		# 1. k_c检测,滑动稳定性检测
		gamma_E1,gamma_G,gamma_E2 = stand["K_c"][type]["gamma_E1"],stand["K_c"][type]["gamma_G"],stand["K_c"][type]["gamma_E2"]
		result_k_c = self.check_k_c(data,gamma_E1,gamma_G,gamma_E2)
		# 2. K_0检测,倾覆稳定性检测
		gamma_E1,gamma_G,gamma_E2 = stand["K_0"][type]["gamma_E1"],stand["K_0"][type]["gamma_G"],stand["K_0"][type]["gamma_E2"]
		result_k_0 = self.check_k_0(data,gamma_E1,gamma_G,gamma_E2)
		# 3. 计算偏心距并检测 e
		e_type = "normal" if type == "I"  else "earthquake"
		e_limit  = stand["e"][self.base_type][e_type]
		result_e = self.check_e(data,e_limit)
		# # 4. sigma检测，基地应力检测
		if type == "I":
			gamma_toe,gamma_heel,gamma_a = stand["sigma"][type]["gamma_toe"],stand["sigma"][type]["gamma_heel"],stand["sigma"][type]["gamma_a"]
		else:
			gamma_toe,gamma_heel,gamma_a = stand["sigma"][type][self.base_type]["gamma_toe"],stand["sigma"][type][self.base_type]["gamma_heel"],stand["sigma"][type][self.base_type]["gamma_a"]
		result_sigma = self.check_sigma(self.parms["base_sigma"],gamma_toe,gamma_heel,gamma_a)
		return {
			"K_c":result_k_c,
			"K_0":result_k_0,
			"e":result_e,
			"sigma":result_sigma
		}
	def get_result(self,img=True):
		# 获取检测结果
		self.cal_exforce()
		#组合I 计算 永久荷载+主可变荷载
		result_check_1 = self.check('I')
		#组合IV 计算 永久荷载+主可变荷载 + 地整荷载
		result_check_4 = self.check('IV')
		return {
			"png_E_3d":self.png_E_3d if img else "",
			"png_stress":self.stress_path,
			"B":self.B, #基底宽度
			"theta":self.theta, #第一破裂角度
			"alpha":self.alpha, #第二破裂角度
			"jx_start":self.to_round(self.jx_start),#第二破裂面
			"jx_end":self.to_round(self.jx_end),#第一破裂面
			"have_2_plm": self.have_2_plm,#是否有第二破裂面
			"area_G":self.to_round(self.area_G),#挡墙截面积
			"area_W":self.to_round(self.area_W),#真实墙背与第二破裂面间土体截面积
			"area_W1":self.to_round(self.area_W1), #第一破裂面与第二破裂面间土体截面积
			"sum_M_y":self.to_round(self.sum_M_y),#总稳定弯矩
			"sum_M_0":self.to_round(self.sum_M_0),#总倾覆弯矩
			"sum_N":self.to_round(self.sum_N),#总竖向力
			"E_a":self.E_a.get_info(), #主动土压力
			"G":self.G.get_info(), #挡墙重力
			"W_toe":self.W_toe.get_info(), #墙趾上土体重力
			"W":self.W.get_info(), #真实墙背与第二破裂面间土体重力
			"W1":self.W1.get_info(), #第一破裂面与第二破裂面间土体重力
			"G_p":self.G_p.get_info(), #挡墙地震力
			"W_p":self.W_p.get_info(), #真实墙背与第二破裂面间土体地震力
			"W1_p":self.W1_p.get_info(), #第一破裂面与第二破裂面间土体地震力
			"check_1":result_check_1, #组合I
			"check_4":result_check_4 #组合IV
		}


####################################################################################
if __name__ == '__main__':
	# 服务器传来的数据
	# 支挡结构设计手册 p179
	# data = {"design":{"H":10,"H2":1,"H3":1,"B":0.5,"B2":0.68,"B3":4.5,"Hzd":1,"L":5.55,"B_fb":0.49,"num_fb":2,"have_tenon":False,"HT":0,"BT":0,"BT1":0,"coords":[[0,0],[0,1],[0.68,1],[0.68,10],[1.18,10],[1.18,1],[5.68,1],[5.68,0],[0,0]],"coord_qz":[5.68,0],"coord_qd":[1.18,10]},"tu":{"po":[[1.18,10],[1.93,10.5],[4.53,10.5],[9.530000000000001,10.5],[14.530000000000001,10.5],[29.53,0.5]],"tuzhu":[{"start":[4.53,10.5],"end":[7.93,10.5],"height":2.7,"width":3.4,"coords":[[4.53,10.5],[4.53,13.2],[7.93,13.2],[7.93,10.5],[4.53,10.5]]},{"start":[9.530000000000001,10.5],"end":[12.930000000000001,10.5],"height":2.7,"width":3.4,"coords":[[9.530000000000001,10.5],[9.530000000000001,13.2],[12.930000000000001,13.2],[12.930000000000001,10.5],[9.530000000000001,10.5]]}]},"parms":{"accuracy":0.5,"f":0.4,"gamma_wall":25,"delta":35,"base_type":"a","gamma_base":19,"base_sigma":400,"phi":35,"gamma_soil":19,"C_z":0.25,"k_H":0.2}}
	# 理正测试数据
	data = {
    "design": {
        "H": 8.5,
        "H2": 0.4,
        "H3": 0.5,
        "B": 0.4,
        "B2": 1.5,
        "B3": 2.5,
        "Hzd": 0.5,
        "L": 6,
        "B_fb": 0.4,
        "num_fb": 2,
        "have_tenon": False,
        "HT": 0,
        "BT": 0,
        "BT1": 0,
        "cover_h": 1.5,
        "coords": [
            [
                0,
                0
            ],
            [
                0,
                0.4
            ],
            [
                1.5,
                0.5
            ],
            [
                1.5,
                8.5
            ],
            [
                1.9,
                8.5
            ],
            [
                1.9,
                0.5
            ],
            [
                4.4,
                0.5
            ],
            [
                4.4,
                0
            ],
            [
                0,
                0
            ]
        ],
        "coord_qz": [
            4.4,
            0
        ],
        "coord_qd": [
            1.9,
            8.5
        ]
    },
    "tu": {
        "po": [
            [
                1.9,
                8.5
            ],
            [
                4.9,
                10.5
            ],
            [
                9.9,
                10.5
            ]
        ],
        "tuzhu": [
            {
                "start": [
                    5.9,
                    10.5
                ],
                "end": [
                    9.100000000000001,
                    10.5
                ],
                "height": 3.1,
                "width": 3.2,
                "coords": [
                    [
                        5.9,
                        10.5
                    ],
                    [
                        5.9,
                        13.6
                    ],
                    [
                        9.100000000000001,
                        13.6
                    ],
                    [
                        9.100000000000001,
                        10.5
                    ],
                    [
                        5.9,
                        10.5
                    ]
                ]
            }
        ]
    },
    "parms": {
        "accuracy": 5,
        "f": 0.5,
        "gamma_wall": 25,
        "base_type": "a",
        "base_sigma": 500,
        "phi": 35,
        "gamma_soil": 19,
        "C_z": 0.25,
        "k_H": 0.2,
        "gamma_0": 1.1
    }
}
	#本地数据库
	# stand = stand

	tmp = Dangtuqiang(data)
	return_data = tmp.get_result(img=False)
	print(return_data)
	print(tmp.sum_M_y)
	print(tmp.sum_M_0)
	print(tmp.sum_N)
	tmp.W_toe.print_data("W_toe")

