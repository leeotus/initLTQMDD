#include "DDpackage.h"


namespace dd {

    void Package::linearSiftingUp(
        short &pos, // variable index
        Edge in, 
        std::map<unsigned short, unsigned short>& Map, 
        std::vector<std::pair<std::pair<short, short>, uint32_t>> &prevMoves
        ) {
        
        //
        const auto n = static_cast<short>(in.p->v + 1);
        prevMoves.push_back({{-1, 0}, size(in)});

        

        while (pos+1 < n ) {
			// std::clog << "向上尝试";
            exchangeBaseCase(pos+1, in, Map); //先对pos，pos-1执行swap
            auto ex_size = size(in);
			linearInPlace(pos+1, in, Map); //再pos，pos-1执行L.T.
			auto lt_size = size(in);
            // std::clog << ex_size <<","<<lt_size<<" ";
            
            if(ex_size <= lt_size){ //swap效果更好
				//** 抵消lt
				// std::clog << "-";
				linearInPlace(pos+1, in, Map);
                auto lt_lt_size = size(in);
                if(lt_lt_size!=ex_size) {
                    std::clog<< lt_lt_size<<", "<<ex_size;
                    throw std::logic_error("lt-lt 后,非还原状态");
                }
                prevMoves.push_back({{SWAP_MOVE, pos+1}, ex_size});

                //debug
                // std::clog<<"swap better"<<std::endl;
			} else { //lt效果更好
                //debug
                // std::clog<<"lt better"<<std::endl;

                prevMoves.push_back({{LINEAR_TRANSFORM_MOVE, pos+1}, lt_size});
			}

            //std::clog << "↓" << ex_size << " ";
            ++pos; //变量位置上移一位 
        }

    }

    void Package::siftingUp(
        short &pos, // variable index
        Edge in, 
        std::map<unsigned short, unsigned short>& Map, 
        std::vector<std::pair<std::pair<short, short>, uint32_t>> &prevMoves
        ) {
        //
        const auto n = static_cast<short>(in.p->v + 1);
        // const auto n = static_cast<short>(in.p->v);
        //
        prevMoves.push_back({{-1, 0}, size(in)});

        while (pos+1 < n ) {
			// std::clog << "向上尝试";
            exchangeBaseCase(pos+1, in, Map); //先对pos，pos-1执行swap
            auto ex_size = size(in);

            prevMoves.push_back({{SWAP_MOVE, pos+1}, ex_size});

            //std::clog << "↓" << ex_size << " ";
            ++pos; //变量位置上移一位 
        }
    }

    void Package::linearSiftingDown(
        short &pos,  
        Edge in, 
        std::map<unsigned short, unsigned short>& Map, 
        std::vector<std::pair<std::pair<short, short>, uint32_t>> &prevMoves
        ) {
        prevMoves.push_back({{-1, 0}, size(in)});

        // std::clog << "向下尝试";
        
       
        while (pos > 0) {
			//std::clog << "向下尝试";
            exchangeBaseCase(pos, in, Map); //先对pos，pos-1执行swap
            auto ex_size = size(in);
			linearInPlace(pos, in, Map); //再pos，pos-1执行L.T.
			auto lt_size = size(in);
            
            if(ex_size <= lt_size){ //swap效果更好
				//** 抵消lt
				// std::clog << "-";
				linearInPlace(pos, in, Map);
                auto lt_lt_size = size(in);
                if(lt_lt_size!=ex_size) {
                    std::clog<< lt_lt_size<<", "<<ex_size;
                    throw std::logic_error("lt-lt 后,非还原状态");
                }
                prevMoves.push_back({{SWAP_MOVE, pos}, ex_size});

                //debug
                // std::clog<<"swap better"<<std::endl;
			} else { //lt效果更好
                

                prevMoves.push_back({{LINEAR_TRANSFORM_MOVE, pos}, lt_size});

                //debug
                // std::clog<<"lt better"<<std::endl;
			}
            

            --pos; //变量位置下移一位 
        }

    }


    void Package::siftingDown(
        short &pos,  
        Edge in, 
        std::map<unsigned short, unsigned short>& Map, 
        std::vector<std::pair<std::pair<short, short>, uint32_t>> &prevMoves
        ) {
        //判断pos与传入的prevMoves是否冲突
        // if(!prevMoves.empty() && prevMoves.back().index-1 != pos){
        //     std::clog<<pos<<", "<< prevMoves.back().index<<";";
        //     throw std::logic_error("传入位置错误,检查sifting down");
        // }

        // std::clog << "向下尝试";
        
        prevMoves.push_back({{-1, 0}, size(in)});

        //

        while (pos > 0) {
			//std::clog << "向下尝试";
            exchangeBaseCase(pos, in, Map); //先对pos，pos-1执行swap
            auto ex_size = size(in);

            prevMoves.push_back({{SWAP_MOVE, pos}, ex_size});

            --pos; //变量位置下移一位 
        }

    }

    int Package::linearSiftingBackward(
        short &optimalPos, 
        Edge in, 
        std::map<unsigned short, unsigned short>& Map, 
        std::vector<std::pair<std::pair<short, short>, uint32_t>> &moves
        ) {
        //找到最小size
        // std::clog<<"回到最好位置 ";
        if(moves.empty()) {
            //
            throw std::logic_error("moves为空.");
        }

        uint32_t bestSize =  moves.front().second;
        for(const auto &i:moves) {
            bestSize = std::min(bestSize, i.second);
        }

        for(auto it=moves.rbegin(); it!=moves.rend(); ) {
            
            if(it->second == bestSize) {
            // if(it == IT) {
                
                // std::clog<<"已是最佳位置a。";
                if(bestSize != size(in)) {throw std::logic_error("best size error"); }
                moves.clear();
                return 1;
            }  
            if(it->first.first == LINEAR_TRANSFORM_MOVE) {
                linearInPlace(it->first.second, in, Map);
            }
           
            exchangeBaseCase(it->first.second, in, Map);
            
            if(it->first.first == INVERSE_TRANSFORM_MOVE) {
                linearInPlace(it->first.second, in, Map);
            }
            optimalPos = it->first.second;
            ++it;
            // if(it->ddsize != size(in)) {
            //     std::clog<<size(in)<<' ';
            //     std::cout<<"还原后的size和记录上的size不一致"<<std::endl;
            //     throw std::logic_error("还原后的size和记录上的size不一致 ");
            // } 
        }

        optimalPos = -1;
        throw std::logic_error("Backward出问题! ");
        return 0;
    }

    // std::vector<Move> Package::undoMoves(
    //     short &pos, 
    //     Edge in, 
    //     std::map<unsigned short, unsigned short>& Map, 
    //     std::vector<Move> &moves
    //     ) {
    //     // std::clog << "回到初始位置: ";

    //     Move invmove;
    //     std::vector<Move> invmoves;
    //     unsigned int curSize = 0;
    //     //
    //     for(std::vector<Move>::reverse_iterator it=moves.rbegin(); it!=moves.rend(); ) {            
    //         if(it->optype == -1) {
    //             if(curSize != it->ddsize) throw std::logic_error("undo() 出错");
    //             break;
    //         }
    //         if(it->optype == SWAP_MOVE) {
    //             //std::clog << "ex" << it->index <<", ";
    //             invmove.optype = SWAP_MOVE;
    //             exchangeBaseCase(it->index, in, Map);
    //         } else if(it->optype == LINEAR_TRANSFORM_MOVE) {
    //             //std::clog << "lt-ex" << it->index <<", ";
    //             invmove.optype = INVERSE_TRANSFORM_MOVE;
    //             linearInPlace(it->index, in, Map);
    //             exchangeBaseCase(it->index, in, Map);
    //         } else {
    //             //std::clog << "ex-lt" << it->index <<", ";
    //             invmove.optype = LINEAR_TRANSFORM_MOVE;
    //             exchangeBaseCase(it->index, in, Map);
    //             linearInPlace(it->index, in, Map);
    //         }
    //         pos = it->pos;
    //         curSize = size(in);
    //         invmove.index = it->index;
    //         invmove.ddsize = curSize;
    //         ++it;
    //         if(it!=moves.rend() && curSize != it->ddsize) {
    //             std::clog<<"recS-curS: "<<it->ddsize<<'-'<<curSize<<' ';
    //             std::cout<<"undo后的size和记录上的size不一致"<<std::endl;
    //             throw std::logic_error("undo后的size和记录上的size不一致 ");
    //         }
            
    //         invmoves.push_back(invmove);
    //     }
    //     return invmoves;
    // }

    Edge Package::dynLinearSifting(Edge in, std::map<unsigned short, unsigned short>& varMap, bool onlySift) {
        const auto n = static_cast<short>(in.p->v + 1); //变量个数
        // const auto n = static_cast<short>(in.p->v);
        

		std::vector<bool> free(n, true); //记录变量是否已经被处理

		computeMatrixProperties = Disabled;
		Edge root{in}; //声明root边


        short pos = -1; //选中变量的index, 电路qubit

		
        static bool matInit = true;
		if(!onlySift && matInit) {
            xorInit(varMap);
            matInit = false;
        }
		

        std::vector<std::pair<std::pair<short, short>, uint32_t>> curMoves;
		for (int i = 0; i < n; ++i) { //遍历各个变量
            assert(is_globally_consistent_dd(in));
            unsigned long min = size(in);
            unsigned long max = 0;

            // std::clog << "    " << i << "/" << n << " size=" << min << " | ";
            for (short j = 0; j < n; j++) {
                if (free.at(varMap[j]) && active.at(varMap[j]) > (unsigned short)max) { //该变量没有被处理过并该变量存在结点
                    max = active.at(varMap[j]); //更新max
                    pos = j; //更新pos，电路变量index
                    assert(max <= std::numeric_limits<int>::max());
                }
            } //到此找到拥有最大结点数的变量，和该变量的索引（位置）
            free.at(varMap[pos]) = false; //设置选中的DD变量为处理状态
            // std::clog<<pos<<','<<varMap[pos]<<" ";
            //move.ddVar = varMap[pos]; //把已经处理的DD变量记录下来
            short optimalPos = pos; 
            short originalPos = pos;

            
            if(pos==n-1) { //选中的变量索引就是顶部
                //
                // std::clog<<"全部向下; "<<std::endl;
                
                // std::clog<<"下: ";
                if(onlySift) siftingDown(pos, in, varMap, curMoves);
                else linearSiftingDown(pos, in, varMap, curMoves);
                // std::clog<<"好: ";
                linearSiftingBackward(optimalPos, in, varMap, curMoves);
                // --optimalPos;
                // std::clog<<"    原来位置和最佳位置："<<originalPos<<"and"<<optimalPos<<"    ";
                
            } else if(pos==0) {
                //
                // std::clog<<"全部向上; "<<std::endl;
                // std::clog<<"上: ";
                if(onlySift) siftingUp(pos, in, varMap, curMoves);
                else linearSiftingUp(pos, in, varMap, curMoves);
                // std::clog<<"好: ";
                linearSiftingBackward(optimalPos, in, varMap, curMoves);
                // std::clog<<"    原来位置和最佳位置："<<originalPos<<"and"<<optimalPos<<"    ";
            } else if(pos < n/2) { // variable is in lower half -> sifting to bottom first
                //
                // std::clog<<"先下后上; "<<std::endl;
                // std::clog<<pos<<' ';
                // std::clog<<"下: ";
                if(onlySift) siftingDown(pos, in, varMap, curMoves);
                else linearSiftingDown(pos, in, varMap, curMoves);
 
                // std::clog<<pos<<' ';
                // std::clog<<"还: ";
                linearSiftingBackward(optimalPos, in, varMap, curMoves);
                if(pos != 0) throw std::logic_error("backwark");
                // std::clog<<pos<<' ';
                // std::clog<<"上: ";
                if(onlySift) siftingUp(pos, in, varMap, curMoves);
                else linearSiftingUp(pos, in, varMap, curMoves);
                // std::clog<<pos<<' ';
                // std::clog<<"好: ";
                linearSiftingBackward(optimalPos, in, varMap, curMoves);
                //--optimalPos;
                // std::clog<<"    原来位置和最佳位置："<<originalPos<<"and"<<optimalPos<<"    ";
               
            } else {
                //
                // std::clog<<"先上后下; "<<std::endl;
                
                // std::clog<<pos<<' ';
                // std::clog<<"上: ";
                if(onlySift) siftingUp(pos, in, varMap, curMoves);
                else linearSiftingUp(pos, in, varMap, curMoves);
                

                linearSiftingBackward(optimalPos, in, varMap, curMoves);
                if(pos != n-1) throw std::logic_error("backwark");
                // --pos; //补偿
                
                // std::clog<<pos<<' ';
                // std::clog<<"下: ";
                if(onlySift) siftingDown(pos, in, varMap, curMoves);
                else linearSiftingDown(pos, in, varMap, curMoves);
                
                // std::clog<<pos<<' ';
                // std::clog<<"好: ";
                linearSiftingBackward(optimalPos, in, varMap, curMoves);
                // --optimalPos;
                // std::clog<<"    原来位置和最佳位置："<<originalPos<<"and"<<optimalPos<<"    ";
              
            }
            
            // if(optimalPos == originalPos) {
            //     // std::clog<<"originalPos-size"<<originalPos<<'-'<<min<<';';
            //     if(min != size(in)) throw std::logic_error("处理后的初始位置的size不等于原初始size");
            // } else if(min < size(in)) {
            //     std::clog << size(in) << ' ';
            //     throw std::logic_error("筛选后结果更差? ");
            // }

            // opSeq.push_back(0);
			initComputeTable();

			// there are nodes which need to renormalized
            if (unnormalizedNodes > 0) {
                // std::clog << "{" << unnormalizedNodes << "} ";
                auto oldroot = root;
                root = renormalize(root);
                decRef(oldroot);
                incRef(root);
                in.p = root.p;
                in.w = root.w;
                if (unnormalizedNodes > 0) {
                    std::clog << "{" << unnormalizedNodes << "} ";
                    throw std::runtime_error("Renormalization failed. " + std::to_string(unnormalizedNodes) + " unnormalized nodes remaining.");
                }
            }
            computeMatrixProperties = Enabled;
            markForMatrixPropertyRecomputation(root); //标记
            recomputeMatrixProperties(root);
            
            //for debug
            // auto scheck = size(in);
            // std::clog << "    " << i << "/" << n << " size=" << scheck << " | ";
		}
		return in; //返回DD指针
    }

    Edge Package::dynSifting(Edge in, std::map<unsigned short, unsigned short>& varMap) {
        const auto n = static_cast<short>(in.p->v + 1); //变量个数
        // const auto n = static_cast<short>(in.p->v);
        // std::clog<<varMap[in.p->v]<<','<<active.at(varMap[in.p->v])<<std::endl;
        

		std::vector<bool> free(n, true); //记录变量是否已经被处理

		computeMatrixProperties = Disabled;
		Edge root{in}; //声明root边


        short pos = -1; //选中变量的index, 电路qubit

		

        std::vector<std::pair<std::pair<short, short>, uint32_t>> curMoveV;
		for (int i = 0; i < n; ++i) { //遍历各个变量
            assert(is_globally_consistent_dd(in));
            unsigned long min = size(in);
            unsigned long max = 0;


            // std::clog << "    " << i << "/" << n << " size=" << min << " | ";
            for (short j = 0; j < n; j++) {
                if (free.at(varMap[j]) && active.at(varMap[j]) > (unsigned short)max) { //该变量没有被处理过并该变量存在结点
                    max = active.at(varMap[j]); //更新max
                    pos = j; //更新pos，电路变量index
                    assert(max <= std::numeric_limits<int>::max());
                }
            } //到此找到拥有最大结点数的变量，和该变量的索引（位置）
            free.at(varMap[pos]) = false; //设置选中的DD变量为处理状态
            // std::clog<<pos<<','<<varMap[pos]<<" ";
            //move.ddVar = varMap[pos]; //把已经处理的DD变量记录下来
            short optimalPos = pos; 
            short originalPos = pos;

            
            if(pos==n-1) { //选中的变量索引就是顶部
                //
                // std::clog<<"全部向下; "<<std::endl;
                
                // std::clog<<"下: ";
                siftingDown(pos, in, varMap, curMoveV);
                if(pos != 0) throw std::logic_error("检查sifting down");
                // std::clog<<"好: ";
                linearSiftingBackward(optimalPos, in, varMap, curMoveV);
                // --optimalPos;
                // std::clog<<"    原来位置和最佳位置："<<originalPos<<"and"<<optimalPos<<"    ";
                
            } else if(pos==0) {
                //
                // std::clog<<"全部向上; "<<std::endl;
                // std::clog<<"上: ";
                siftingUp(pos, in, varMap, curMoveV);
                if(pos != n-1) throw std::logic_error("检查sifting up");
                // std::clog<<"好: ";
                linearSiftingBackward(optimalPos, in, varMap, curMoveV);
                // std::clog<<"    原来位置和最佳位置："<<originalPos<<"and"<<optimalPos<<"    ";
            } else if(pos < n/2) { // variable is in lower half -> sifting to bottom first
                //
                // std::clog<<"先下后上; "<<std::endl;
                // std::clog<<pos<<' ';
                // std::clog<<"下: ";
                siftingDown(pos, in, varMap, curMoveV);
                if(pos != 0) throw std::logic_error("检查sifting down 2");
                siftingUp(pos, in, varMap, curMoveV);
                // std::clog<<pos<<' ';
                // std::clog<<"好: ";
                linearSiftingBackward(optimalPos, in, varMap, curMoveV);
                //--optimalPos;
                // std::clog<<"    原来位置和最佳位置："<<originalPos<<"and"<<optimalPos<<"    ";
               
            } else {
                //
                // std::clog<<"先上后下; "<<std::endl;
                
                // std::clog<<pos<<' ';
                // std::clog<<"上: ";
                siftingUp(pos, in, varMap, curMoveV);
                if(pos != n-1) throw std::logic_error("检查sifting up");
                
                // --pos; //补偿
                
                // std::clog<<pos<<' ';
                // std::clog<<"下: ";
                siftingDown(pos, in, varMap, curMoveV);
                
                // std::clog<<pos<<' ';
                // std::clog<<"好: ";
                linearSiftingBackward(optimalPos, in, varMap, curMoveV);
                // --optimalPos;
                // std::clog<<"    原来位置和最佳位置："<<originalPos<<"and"<<optimalPos<<"    ";
              
            }
            
            // if(optimalPos == originalPos) {
            //     // std::clog<<"originalPos-size"<<originalPos<<'-'<<min<<';';
            //     if(min != size(in)) throw std::logic_error("处理后的初始位置的size不等于原初始size");
            // } else if(min < size(in)) {
            //     std::clog << size(in) << ' ';
            //     throw std::logic_error("筛选后结果更差? ");
            // }

            // opSeq.push_back(0);
			initComputeTable();

			// there are nodes which need to renormalized
            // if (unnormalizedNodes > 0) {
            //     // std::clog << "{" << unnormalizedNodes << "} ";
            //     auto oldroot = root;
            //     root = renormalize(root);
                
            //     decRef(oldroot);    
            //     incRef(root);
            //     in.p = root.p;
            //     in.w = root.w;
                
            //     // std::clog << "{" << "renormalize" << "} ";
                
            //     if (unnormalizedNodes > 0) {
            //         std::clog << "{" << unnormalizedNodes << "} ";
            //         // throw std::runtime_error("Renormalization failed. " + std::to_string(unnormalizedNodes) + " unnormalized nodes remaining.");
            //     }
            // }
            computeMatrixProperties = Enabled;
            markForMatrixPropertyRecomputation(root); //标记
            recomputeMatrixProperties(root);
            
            //for debug
            // auto scheck = size(in);
            // std::clog << "    " << i << "/" << n << " size=" << scheck << " | ";
		}
		return in; //返回DD指针
    }

    // Edge Package::dynFlatLinearSifting(Edge in, std::map<unsigned short, unsigned short>& varMap, std::deque<Edge> &allDD, bool onlySift) {
    //     const auto n = static_cast<short>(in.p->v + 1); //变量个数
    //     //Move move;
    //     std::vector<Move> moveUp; //list of up moves
    //     std::vector<Move> moveDown;

	// 	std::vector<bool> free(n, true); //记录变量是否已经被处理

	// 	computeMatrixProperties = Disabled;
	// 	Edge root{in}; //声明root边


    //     short pos = -1; //选中变量的index, 电路qubit

		
    //     static bool matInit = true;
	// 	if(!onlySift && matInit) {
    //         xorInit(varMap);
    //         matInit = false;
    //     }
		


	// 	for (int i = 0; i < n; ++i) { //遍历各个变量
    //         assert(is_globally_consistent_dd(in));
    //         unsigned long min = size(in);
    //         unsigned long max = 0;

    //         // std::clog << "    " << i << "/" << n << " size=" << min << " | ";
    //         for (short j = 0; j < n; j++) {
    //             if (free.at(varMap[j]) && active.at(varMap[j]) > (unsigned short)max) { //该变量没有被处理过并该变量存在结点
    //                 max = active.at(varMap[j]); //更新max
    //                 pos = j; //更新pos，电路变量index
    //                 assert(max <= std::numeric_limits<int>::max());
    //             }
    //         } //到此找到拥有最大结点数的变量，和该变量的索引（位置）
    //         free.at(varMap[pos]) = false; //设置选中的DD变量为处理状态
    //         //move.ddVar = varMap[pos]; //把已经处理的DD变量记录下来
    //         short optimalPos = pos; 
    //         short originalPos = pos;

    //         //记录初始位置
    //         Move curState;
    //         curState.ddsize = min;
    //         curState.index = -1;
    //         curState.pos = pos;
    //         curState.optype = -1;
    //         std::vector<Move> curMoveV;
    //         curMoveV.push_back(curState);
            
    //         if(pos==n-1) { //选中的变量索引就是顶部
    //             //
    //             // std::clog<<"全部向下; "<<std::endl;
                
    //             // std::clog<<"下: ";
    //             if(onlySift) siftingDown(pos, in, varMap, curMoveV);
    //             else linearSiftingDown(pos, in, varMap, curMoveV);
    //             // std::clog<<"好: ";
    //             linearSiftingBackward(optimalPos, in, varMap, moveDown);
    //             // --optimalPos;
    //             // std::clog<<"    原来位置和最佳位置："<<originalPos<<"and"<<optimalPos<<"    ";
                
    //         } else if(pos==0) {
    //             //
    //             // std::clog<<"全部向上; "<<std::endl;
    //             // std::clog<<"上: ";
    //             if(onlySift) siftingUp(pos, in, varMap, curMoveV);
    //             else linearSiftingUp(pos, in, varMap, curMoveV);
    //             // std::clog<<"好: ";
    //             linearSiftingBackward(optimalPos, in, varMap, moveUp);
    //             // std::clog<<"    原来位置和最佳位置："<<originalPos<<"and"<<optimalPos<<"    ";
    //         } else if(pos < n/2) { // variable is in lower half -> sifting to bottom first
    //             //
    //             // std::clog<<"先下后上; "<<std::endl;
    //             // std::clog<<pos<<' ';
    //             // std::clog<<"下: ";
    //             if(onlySift) siftingDown(pos, in, varMap, curMoveV);
    //             else linearSiftingDown(pos, in, varMap, curMoveV);
    //             if(moveDown.back().index != pos+1) {
    //                 throw std::logic_error("检查aux-lsDown");
    //             }
    //             // std::clog<<pos<<' ';
    //             // std::clog<<"还: ";
    //             undoMoves(pos, in, varMap, moveDown);
    //             if(min != size(in)) throw std::logic_error("undo 错误 ");
    //             // std::clog<<pos<<' ';
    //             // std::clog<<"上: ";
    //             if(onlySift) siftingUp(pos, in, varMap, moveUp);
    //             else linearSiftingUp(pos, in, varMap, moveUp);
    //             // std::clog<<pos<<' ';
    //             // std::clog<<"好: ";
    //             linearSiftingBackward(optimalPos, in, varMap, moveUp);
    //             //--optimalPos;
    //             // std::clog<<"    原来位置和最佳位置："<<originalPos<<"and"<<optimalPos<<"    ";
               
    //         } else {
    //             //
    //             // std::clog<<"先上后下; "<<std::endl;
                
    //             // std::clog<<pos<<' ';
    //             // std::clog<<"上: ";
    //             if(onlySift) siftingUp(pos, in, varMap, curMoveV);
    //             else linearSiftingUp(pos, in, varMap, curMoveV);
    //             // if(moveUp.back().index != pos) {
    //             //     throw std::logic_error("检查aux-lsUp");
    //             // }
                
    //             // std::clog<<pos<<' ';
    //             // std::clog<<"还: ";
    //             undoMoves(pos, in, varMap, moveUp);
    //             if(min != size(in)) throw std::logic_error("undo 错误 ");
    //             // --pos; //补偿
                
    //             // std::clog<<pos<<' ';
    //             // std::clog<<"下: ";
    //             if(onlySift) siftingDown(pos, in, varMap, moveDown);
    //             else linearSiftingDown(pos, in, varMap, moveDown);
                
    //             // std::clog<<pos<<' ';
    //             // std::clog<<"好: ";
    //             linearSiftingBackward(optimalPos, in, varMap, moveDown);
    //             // --optimalPos;
    //             // std::clog<<"    原来位置和最佳位置："<<originalPos<<"and"<<optimalPos<<"    ";
              
    //         }
            
    //         if(optimalPos == originalPos) {
    //             // std::clog<<"originalPos-size"<<originalPos<<'-'<<min<<';';
    //             if(min != size(in)) throw std::logic_error("处理后的初始位置的size不等于原初始size");
    //         } else if(min < size(in)) {
    //             std::clog << size(in) << ' ';
    //             throw std::logic_error("筛选后结果更差? ");
    //         }

    //         // opSeq.push_back(0);
	// 		initComputeTable();

	// 		// there are nodes which need to renormalized
    //         if (unnormalizedNodes > 0) {
    //             std::clog << "{" << unnormalizedNodes << "} ";
    //             auto oldroot = root;
    //             root = renormalize(root);
    //             incRef(root);
    //             decRef(oldroot);
    //             in.p = root.p;
    //             in.w = root.w;
    //             std::clog << "{" << unnormalizedNodes << "} ";
    //             if (unnormalizedNodes > 0) {
    //                 uint32_t deqNum = allDD.size();
    //                 for(; deqNum>0; --deqNum) {
    //                     root = allDD.front();
    //                     allDD.pop_front();
    //                     auto oldroot = root;
    //                     root = renormalize(root);
    //                     decRef(oldroot);
    //                     incRef(root);
    //                     allDD.push_back(root);
    //                 }
    //                 std::clog << "{" << unnormalizedNodes << "} ";
    //                 throw std::runtime_error("Renormalization failed. " + std::to_string(unnormalizedNodes) + " unnormalized nodes remaining.");
    //             }
    //         }
    //         computeMatrixProperties = Enabled;
    //         markForMatrixPropertyRecomputation(root); //标记
    //         recomputeMatrixProperties(root);
            
    //         //for debug
    //         // auto scheck = size(in);
    //         // std::clog << "    " << i << "/" << n << " size=" << scheck << " | ";
	// 	}
	// 	return in; //返回DD指针
    // }



}