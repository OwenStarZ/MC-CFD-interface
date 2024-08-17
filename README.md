# 接口说明与使用规范
---
## 接口说明
本接口程序能够进行蒙卡程序RMC与CFD程序ANSYS Fluent的稳态、瞬态核热耦合计算，目前只适用于Linux平台，程序中每个文件对应不同的功能如下：
* \lib_h5rw文件夹为外部挂载hdf5读写所需动态链接库的编译路径
* \libudf文件夹为Fluent加载的udf编译路径与库文件，Windows平台下为.dll与.lib；Linux平台下为.so
* compile.sh为在Linux平台下重编译上述两个动态链接库的脚本
* init.py为耦合初始化脚本
* iter.py为耦合迭代控制脚本
* runFluent.jou为Fluent控制脚本，文件名中最后一位数字为0的为稳态控制脚本，1为瞬态控制脚本
* **coupling.dat为【耦合输入卡】**
  #### 耦合输入卡说明
  * 第1行为稳态或瞬态单时间步内的**最大Picard迭代次数**
  * 第2行为稳态或**当前**瞬态时间步内的模型**总功率**，单位为MW
  * 第3行为光子等直接产热份额，一般设为0.0即可
  * 第4~6行为Fluent几何模型相对于RMC几何模型的笛卡尔坐标xyz偏移值，单位为cm，建议建模时就做到几何对应，则这几行数据均为0.0
  * 第7~10行为默认，也是稳态情形下初始的燃料(fuel)，慢化剂(moderator)，反射层(reflector)，冷却剂(coolant)的温度，第11行为默认冷却剂密度，均为SI单位制
  * #Project Name为RMC输入卡与Fluent case的文件名（不含后缀）
  * #Coupling mode为计算模式，0代表稳态，1代表瞬态。*#大时间步瞬态，即燃耗计算尚在开发中*
  * #time step为时间步设置，**非常重要**，稳态情形下只有第五个量起作用，瞬态情形则都有作用：
    * 第一个量为**瞬态核热耦合**的时间步，又称外迭代时间步，目前只支持固定时间步
    * 第二个量为**瞬态核热耦合**的总时长
    * 第三个量为**中子物理场**的点堆参数更新子时间步，又称中子内迭代时间步
    * 第四个量为**热工水力场**的瞬态方程时间步，又称热工水力内迭代时间步
    * 第五个量为热工水力场稳态方程或瞬态方程每个时间步的最大迭代次数
  * #DIM到#OUTDIMR分别为功率网格，燃料、冷却剂、慢化剂、反射层热工参数信息传输网格分别在笛卡尔坐标xyz方向上的划分数量，目前只支持均匀结构化划分
  * #P_bndry为功率网格的边界位置，**单位为cm**
  * #F_bndry ~ #R_bndry分别为不同材料热工参数信息传输网格的边界位置，**单位为m**
  * #Multilevel Flag为用01控制的多尺度模型开关，常应用于含有TRISO颗粒的模型，如HTGR
  * #Omega为耦合传输松弛因子，常取0.5，**目前默认功率一定使用松弛法**，如不想松弛，需指定该值为1
  * #SOR flag for th_calculation为是否对热工数据使用松弛法，如使用，与功率一致。**目前的研究表明使用的意义不大，建议关闭**
  * #Residual factor为功率的上α分位点，通常取3.0即可，如模型复杂可再度适当放宽
  * #tmax discrepancy tolerance为最大温度相对误差的容许上限
  * #Parallel processes为RMC与Fluent的并行**进程数**这两个数字建议均不要超过计算机的**核数，而非线程数**
  * #Output all iter flag为是否无条件输出所有迭代步，即无需判敛
  * #Coupling order为耦合顺序，**只对瞬态有效**，0为先从RMC开始，1为先从Fluent开始，默认为0

使用接口程序之前，需要进行环境变量配置，输入命令打开环境变量配置：
```
vim ~/.bashrc
```
单击键盘``i``键进入编辑模式，输入以下内容并替换中括号为实际数字或路径：
```
ulimit -s [number]                   # 取消栈内存限制，[number]数字单位为kbytes
export PATH=[PATH to mpich]/mpich/bin:$PATH
export PATH=[PATH to hdf5]/hdf5/[version number]:$PATH
export PATH=[PATH to ansys_inc]/ansys_inc/v211/fluent/bin:$PATH
export PATH=[PATH to RMC folder]/RMC:$PATH
export PATH=[PATH to mpich]/ mpich/bin:$PATH
export LD_LIBRARY_PATH=[PATH to hdf5]/hdf5/[version number]/lib:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=[PATH to hdf5]/hdf5/[version number]/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=[PATH to hdf5]/hdf5/[version number]/include:$CPLUS_INCLUDE_PATH
export RMC_DATA_PATH=[PATH to RMC folder]/RMC
```
最后输入命令以应用配置：
```
source .bashrc
```

---

## 使用规范

使用时先将本接口程序下的所有文件复制到工作目录下。**此外还需复制RMC截面库索引文件xsdir到工作目录**，本接口程序不进行提供

### 稳态情形

* **解决单物理场**。准备好RMC输入卡与**不挂载UDF的**Fluent case文件
    * **中子物理场**：准备.rmc格式的输入卡
      * 中子代信息**必须手动设置**，所有Tally**必须手动设置**，且MeshTally1必须将网格信息单独起一行，其余Tally禁止与MeshTally1采用相同换行策略
      * 使用Plot功能检查中子物理模型几何的正确性，注意要修改需要考虑耦合效应的栅元的**几何重复策略**
      * 先手动指定固定的栅元温度，检查是否能正常进行临界计算
      * 修改coupling.dat输入卡中热工参数传输网格部分的设置，在工作目录打开命令行，输入
        ```
        python init.py
        ```
        生成.h5格式的试验热工参数分布文件，同时按RMC使用规范修改输入卡，检查能否正常读取外部热工参数文件并进行临界计算，这一步的目的是检查模型是否具备耦合条件，以及占用内存是否溢出。可以输入
        ```
        watch free -g
        ```
        监控内存占用情况
      * 全部通过后，重命名.rmc输入卡，名称（不含后缀）需与coupling.dat中的项目名称对应
  * **热工水力场**：准备.cas.h5格式的算例文件
    * 先准备好几何模型与网格模型的工程文件
    * 导入Fluent求解器，先不挂载任何UDF，使用**恒定热源**，设置好边界等条件，检查是否能正常计算，以检查模型正确性
    * 【可选】如果有UDF定义的物性，编译并挂载物性相关UDF，检查是否能正常计算，以检查物性UDF编译正确性
    * 做网格无关性和湍流壁面函数*y*^+^的验证，对Standard *k*-$\varepsilon$ model，应满足30<*y*^+^<200，对Standard *k*-$\omega$ model，应满足*y*^+^<5
    * 全部通过后保存case，**这个不挂载UDF（不包含物性UDF）的case建议单独保存备用**
* GUI模式下打开.cas.h5，记录所有要考虑耦合效应的，不同材料对应区域的**zone ID**，修改udf.c中的对应行【这一步非常非常重要！】
* 打开coupling.dat输入卡，**逐行检查耦合计算设置并修改**
* 打开iter.py，**仅需要**修改需要归档的文件名称
* 目前，稳态耦合工作目录下应该包含该项目的所有接口文件，以及.rmc与.cas.h5文件，此时输入
  ```
  sh compile.sh
  ```
  会生成info_fuel.h5等四个初始热工参数文件，以及动态链接库libh5rw.so
* **手动挂载**libudf，**另存**.cas.h5并退出，**另存后的文件名称需与耦合输入卡中的名称一致**
* 输入
  ```
  python iter.py
  ```
---
### 瞬态情形

* **先进行稳态耦合计算**（建议在其它目录下进行避免混淆）。能够进行稳态耦合计算说明模型几何，文件读写，内存占用，以及UDF都不存在问题。能够进行瞬态计算的**前提是稳态耦合的*k*~eff~收敛值与1相差1~2个标准差之内**
* 将稳态耦合中的接口文件复制到瞬态耦合工作目录下
* **解决单物理场**。准备好RMC与Fluent的输入文件
  * **中子物理场**：准备.rmc，.innerproduct，.State.h5格式的输入卡
    * 满足*k*~eff~的条件时，在稳态耦合的.rmc输入卡的基础上修改，将临界计算改为稳态时空动力学计算，即QUASISTATIC_S计算，使用收敛后的热工参数文件。**修改外迭代时间步长与中子数**，得到.innerproduct，.State.h5格式的文件
    * 将QUASISTATIC_S输入卡复制到瞬态工作目录下，将QUASISTATIC_S修改为QUASISTATIC_D，并加入内迭代步数信息（此处内迭代信息可以任意指定，脚本会根据coupling.dat中的信息修改这一数值），**但禁中子数与外迭代步长必须与QUASISTATIC_S的设置一致**。同时，也将.innerproduct，.State.h5格式的文件复制到瞬态耦合工作目录，由于.State.h5会在瞬态计算流程中被反复覆盖，**须将.State.h5格式文件重命名**，例如Singlepin_init.State.h5
    * 打开coupling.dat，按照外迭代时间步长度，修改QUASISTATIC_D的内迭代时间步数量，使得内迭代时间步长度保持不变，例如，对外迭代时间步长分别为0.5s或1s的算例，则分别修改内迭代时间步数量为50与100，则内迭代时间步长均为0.01s。*#耦合外迭代与中子内迭代时间步长的选取尚待研究*
  * **热工水力场**：准备.cas.h5格式的算例文件，以及.dat.h5格式的初始条件
    * 在GUI中将**不挂载任何UDF**的稳态算例导入Fluent求解器，并导入收敛后的稳态CFD结果.dat.h5，使用**恒定热源，并将稳态改为瞬态模式**，检查是否能正常计算，以检查模型正确性
    * 对内迭代时间步长做Courant数的验证，根据ANSYS Fluent user manual的说明，最大的Courant数不应超过20~40这一范围，且每个内迭代时间步应在10次内收敛。理想条件下，Courant数越接近1越好
    * 全部通过后，**编译挂载UDF，并保存**。由于瞬态计算在稳态计算的基础上进行，所以不应修改任何数据传输网格的划分策略，因此也无需修改UDF，直接编译挂载即可。需要注意的是，这一个重复挂载的过程建议一定要做，避免UDF路径的混淆
    * 将稳态耦合收敛后的.dat.h5格式的CFD结果文件复制到瞬态耦合工作目录**用作瞬态耦合的初始条件**，由于.dat.h5会在瞬态计算流程中被反复覆盖，**须将.dat.h5格式文件重命名**，例如pin_init.dat.h5
* 将稳态耦合所得收敛的功率，热工参数传输文件，即MeshTally1.h5，以及info_fuel.h5等四个文件复制到瞬态耦合工作目录下
* 打开coupling.dat输入卡，逐行检查耦合计算设置并修改，**特别注意要修改第二行功率值为初始功率**
* 打开iter.py，**仅需要**修改需要归档的文件名称
* 目前，瞬态耦合工作目录下应该包含该项目的所有接口文件，.rmc与.cas.h5文件，以及包含初始条件的八个文件，包括：.innerproduct; _init.State.h5; MeshTally1.h5; _init.dat.h5; info_fuel.h5等四个文件，此时输入
  ```
  sh compile.sh
  ```
  会生成动态链接库libh5rw.so
* 输入
  ```
  python iter.py
  ```
