# matlab-actions
Run matlab scripts in github-actions.

使用GitHub服务器跑耗时的MATLAB代码。

## 1.MATLAB代码修改
在MATLAB主函数代码最后使用命令保存所需数据，如：

    save data.mat;

    print('-f1','-dpng','savepic1.png'); 
## 2.创建工作区
在根目录创建工作区，如：

    workspace1

将MATLAB脚本放入工作区中

## 3.获取GitHub_Token
头像旁边的箭头-->Settings-->Developer settings-->Personal access tokens-->Generate new token

![image](https://user-images.githubusercontent.com/63141816/169046576-673bf3e9-3427-4096-a43e-d23917286188.png)

## 4.将Token填入Actions
Settings-->Secrets-->Actions-->New repository secret

![image](https://user-images.githubusercontent.com/63141816/169047991-cc90d5bd-2c66-45ce-9898-aba1acfc92a9.png)

## 5.开始运行Actions
根据提示更改运行参数

![image](https://user-images.githubusercontent.com/63141816/169052140-4cdcd486-19a2-4e55-9335-2a2f3d206d2a.png)

然后点击

    Run workflow

等待运行结束，数据将会发布在Releases处

## 特别提醒：运行重要代码时，请将仓库设为私密

# 鸣谢
* [MATLAB Actions](https://github.com/matlab-actions)

