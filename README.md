# UAV_Integrated_Navigation




## Preface

1. 遵循“先解释，后执行”的原则，注释总是写在对应程序的前面（上方）

2. 注意：在matlab程序中表示同时具有上标和下标的变量时，严老师书上的程序案例按照“上标在前，下标在后”的原则，我的习惯是按照 “下标在前，上标在后”的原则。eg：书中对于 "从b系到n系的过渡矩阵" Cbn 的写法为Cnb。

3. 程序中,总是按照如下方式定义变量：
   1. pos : [ 纬度, 经度, 高度 ] 单位为rad, rad, m
   2. vn : [ vE, vN, vU ] 单位为m/s
   
   若无特地声明，表示角的大小的变量的单位总是rad，其他物理量的单位总是相应的国际单位。
   
4. 以_ref结尾的变量，注释中称其为“参考值”，“标准值”。这两种说法均表示“真实值”或“没有误差的值”。其中，本人更愿意把实测数据中的高精度参考称为“参考值”，把仿真程序生成的绝对准确的结果称为“标准值”，但真正写到程序里时也没区分那么细。

5. 通过仿真程序生成的数据以_SD作为结尾。

6. k_init 特指初始时刻数据的对应编号



## Meaning of some abbreviation

**ref**: reference    **msr**: measurement    **sv**: save    **err**: error    **prv**: previous    **opt**: optimal   **sml**: simulation    **crt**: current    **init**: initial    **num**: number    **smpl**: sample    **cmpt**: compute    **Lp**: loop    **idx**: index    **SD**: simulation data  

**导航专用缩写:**

**avp**: attitude、velocity、position    **pos**: positon    **att**: attitude    **vn**: ground velocity     **acc**: accelerometer    **gyr**: gyroscope    **ang**: angle    



## What's included in this code Repo.?

- 惯导解算程序
- 卡尔曼间接滤波算法（导航误差作为状态变量）
  - 输出校正：估计的误差状态仅用来更新导航输出，不对惯导解算进行纠正
  - 反馈校正：利用误差估计对惯导解算过程进行校正



## Further work

* [ ] Sage-Husa自适应滤波算法
* [ ] 强跟踪卡尔曼滤波算法
* [ ] ... ...



## References

- [1] 严恭敏, 捷联惯导算法与组合导航原理.
- [2] 王辰熙, 南航
