# Top K Shortest Paths 问题

## Djikstra算法

Dijkstra算法可以解决什么问题呢? —— "单源最短路径"。

也就是说, 给定一个起点, 我们可以利用Dijkstra算法求出"起点到其他点的最短路径"。


## Top K Shortest Paths 简介

现在我想要处理这么一个问题:
```math
\quad 给定起点, \textbf{终点}
```

很自然地, 我们可以利用Dijkstra算法求出"起点到终点的最短路径"。

但是有这么一个问题:

```math
\text{这个最短路径只有} \underline{1} \text{条}
```

如果我想要:

```math
\text{起点到终点的} \underline{K} \text{条最短路径}
```

那怎么才能处理呢? —— "Djikstra Algo" 和 "Yen Algo"。

## 限制条件

如果我已经处理好了Top K Shortest Paths, 那我还想要再添加一些限制条件, 就比如说:

1. Ban掉某个点 S.t. 路径完全不经过该点
2. 相应地, 会有一个Unban, 将Ban掉的点从中解除
3. MaxTrans:  限制最大中继点个数
4. paths:  设定输出的最短路径个数

问题在于:

1. 如何保持整个算法的简洁性
2. 如何保持整个算法的内在操作逻辑尽可能地简单
