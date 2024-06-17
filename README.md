# 接口说明与使用规范：

## 接口说明：

* \lib_h5rw文件夹为外挂hdf5读写所需动态链接库的编译路径
* \libudf文件夹为Fluent加载的udf编译路径与库文件
* compile.sh为在Linux平台下重编译上述两个动态链接库的脚本
* coupling.dat为【耦合输入卡】
* init.py为耦合初始化脚本
* iter.py为耦合迭代控制脚本
* runFluent.jou为Fluent控制脚本

## 使用规范：

### 首次使用，或需要修改数据传递网格尺度的情况下：

* 先解决单物理场。准备好RMC输入卡与Fluent case的.rmc与.cas文件【设好RMC中子代信息，Mesh信息，与.jou文件中的迭代步信息】
【此处建议先在外部准备好一个未挂载libudf的case，再复制到工作路径下】
* 打开.cas文件，检查燃料、冷却剂等反馈作用区域对应的zone ID，并修改udf.c中的对应行【这一步对同一几何只需要做一次】
* 打开.rmc输入卡，修改对应位置的重复几何策略，需要与【热工数据统计网格】的尺度尽可能对应，【无需修改其它】
* 打开.dat输入卡，【必须逐行检查耦合计算设置并修改！！！！！！！！】
* 打开iter.py，【仅需要】修改需要归档的文件名称
* sh compile.sh
* 【手动挂载】libudf，保存.cas并退出
* python iter.py

### 无需重新编译libudf情况下（数据传递网格尺度不变，多尺度计算开关状态不变）：

* 先解决单物理场。准备好RMC输入卡与Fluent case的.rmc与.cas文件【设好RMC中子代信息; Mesh信息与.jou文件中的迭代步信息会自动更新】
    * 稳态情形下需要.rmc与.case.h5即可，init.py会提供初始温度分布供RMC计算，而后由首次RMC计算得到首次CFD计算所需的功率分布
    * 瞬态情形下除.rmc与.case.h5以外，尚需_init.rmc.State.h5与_init.dat.h5，
      除此外需要info_fuel.h5等四个初始热工参数文件，与MeshTally.h5、.innerproduct的初始功率形状、幅度文件。
      这些文件应由此前运行过的稳态或上一个大时间步的瞬态计算得到。
* 打开coupling.dat输入卡，【必须逐行检查耦合计算设置并修改！！！！！！！！】【对瞬态情形，不建议修改网格尺度】
* python init.py
* python iter.py