

MyGo <- function(n=19)
{
  ###初始设置
  ##构造棋盘
  if (!interactive()) return()
  par(mar = rep(0, 4))
  board <- function()
  {
    plot(1:n, type = "n", xlim = c(1, n), axes = FALSE, xlab = "",
         ylab = "", bty = "o", lab = c(n, n, 1))
    segments(1, 1:n, n, 1:n)
    segments(1:n, 1, 1:n, n)
    points(rep(c(4, 10, 16), 3), rep(c(4, 10, 16), each = 3),
           pch = 19, cex = 1.2)
  }    
  board()
  box()
  
  ##三大列表初始设置
  playedlist <- c()    #历史记录
  conset <- list(list(),list())    #连通分支列表,是二重嵌套列表，第一子列表属黑棋，第二属白棋
  bouset <- list(list(0+0i),list(0+0i))    #连通分支边界列表，同上
  ncb=c(0,0)    #分支数
  ntake <- c(0,0);eat <- 0
  
  ##棋盘特殊位置
  SW<-1+1i    #棋盘四角
  SE<-19+1i
  NE<-19+19i
  NW<-1+19i
  S<- 2:18+1i    #棋盘四边（不含角）
  E<- complex(0,19,2:18)
  N<- 2:18+19i
  W<- complex(0,1,2:18)
  dire<-c(1i,-1i,-1,1)    #棋盘四合
  
  ##边界函数Bd（从中去掉历史位置即可得到“气”）
  Bd<-function(A)
  {
    B=c()
    for(a in A)
    {
      if(a==SW){
        pa=c(1+2i,2+1i)
      }else if(a==SE){
        pa=c(18+1i,19+2i)
      }else if(a==NE){
        pa=c(18+19i,19+18i)
      }else if(a==NW){
        pa=c(1+18i,2+19i)
      }else if(is.element(a,S)){
        pa=c(a-1,a+1,a+1i)
      }else if(is.element(a,E)){
        pa=c(a-1,a-1i,a+1i)
      }else if(is.element(a,N)){
        pa=c(a-1,a+1,a-1i)
      }else if(is.element(a,W)){
        pa=c(a+1,a-1i,a+1i)
      }else{
        pa <- a+dire
      }
      for(i in pa) if(!is.element(i,A)) B=c(B,i)
    }
    B
  }
  
  ##提子
  take <- function()    #把被提的子从分支列表和历史记录删去后，将棋盘重新绘制
  {
    board()
    BB <- unlist(conset[[1]])
    WW <- unlist(conset[[2]])
    points(Re(BB),Im(BB), cex = 3, pch = 19, bg = "black")
    points(Re(WW),Im(WW), cex = 3, pch = 21, bg = "white")
  }
  
  ###进入对弈
  k <- 0
  repeat    
  {
    for (j in 1:2)    ##黑白交替落子
    {
      repeat
      {
        l <- locator(1)    #获得落子坐标
        l$x <- min(n, max(1, round(l$x)))
        l$y <- min(n, max(1, round(l$y)))
        xy <- complex(0,l$x,l$y)
        if (!is.element(xy, playedlist))    #禁走历史位置
          break
      }
      #落子
      points(l, cex = 3, pch = c(19, 21)[j], bg = c("black", "white")[j])
      print(c(l$x,l$y))
      
      
      ##损气提吃
      s=1
      while(s<=ncb[3-j])
      {  
        opp <- bouset[[3-j]][[s]]    #某方第 s 个分支的边界              
        if(is.element(xy,opp))
        {
          if(length(opp)==1)    #待提的分支们
          {
            playedlist <- setdiff(playedlist,conset[[3-j]][[s]])
            ntake[j] <- ntake[j]+length(conset[[3-j]][[s]]);eat <- 1
            bouset[[3-j]][[s]] <- NULL    #删除被提的分支记录
            conset[[3-j]][[s]]<-NULL
            ncb[3-j] <- ncb[3-j]-1    #对方分支数减少
            for(t in 1:ncb[j])
            {
              bouset[[j]][[t]] <- setdiff(Bd(conset[[j]][[t]]),playedlist)
            }
            s <- s-1
          }else{
            bouset[[3-j]][[s]] <- setdiff(opp,xy)      
          }
        }
        s <- s+1
      }
      
      ###构造三大列表（历史、连通分支、边界）
      p=0    #分支构建指标(p>0分支合并、p=0新建分支)
      U <- c()    #合并分支
      i <- 1;r<- c()
      while(i<=ncb[j]) 
      {
        if(ncb[j]==0) break
        if(is.element(xy, bouset[[j]][[i]]))    #新下的棋属于哪个连通分支
        {
          p <- p+1
          if(p==1)
          {
            q <- i
            U <- conset[[j]][[i]]
            conset[[j]][[i]]<-NULL
            bouset[[j]][[i]]<-NULL
            ncb[[j]] <- ncb[[j]]-1
            i <- i-1
          }
          if(p>1)
          {
            U <- union(U,conset[[j]][[i]])
            conset[[j]][[i]]<-NULL
            bouset[[j]][[i]]<-NULL
            i=i-1
            ncb[j] <- ncb[j]-1 
          }
        }
        i <- i+1
      }
      
      ncb[j] <- ncb[j]+1
      if(p==0){
        if(ncb[j]==0)
        {
          conset[[j]] <- list(xy)    #孤子建新支
          bouset[[j]] <- list(setdiff(Bd(xy),playedlist))
        }else{
          conset[[j]][[ncb[j]]] <- xy    
          bouset[[j]][[ncb[j]]] <- setdiff(Bd(xy),playedlist)
        }
      }else if(p==1){
        if(ncb[j]>0)
        {
          conset[[j]][[ncb[j]]] <- c(U,xy)
          bouset[[j]][[ncb[j]]] <- setdiff(Bd(c(U,xy)),playedlist)
        }else{
          conset[[j]] <- list(c(U,xy))
          bouset[[j]] <- list(setdiff(Bd(c(U,xy)),playedlist))
        }
      }else{
        conset[[j]][[ncb[j]]] <- c(U,xy)
        bouset[[j]][[ncb[j]]] <- setdiff(Bd(c(U,xy)),playedlist)
      }
      
      playedlist <- c(playedlist, xy)
      if(eat>0)
      {
        print(list(黑方提子数=ntake[1],白方提子数=ntake[2]));eat <- 0
        take()
      }
      
      k <- k+1
      if (k >= n^2)break    ###满盘结束(满盘皆输555~)
    }
    if(k >= n^2) break
  }
}
MyGo()
