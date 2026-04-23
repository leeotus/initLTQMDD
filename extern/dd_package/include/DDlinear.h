#ifndef DDlinear_H
#define DDlinear_H
#include <stack>

#define SWAP_MOVE 0
#define LINEAR_TRANSFORM_MOVE 1
#define INVERSE_TRANSFORM_MOVE 2
#define UPDATE_MAP 3
#define LOW_LINEAR_TRANSFORM_MOVE 5
#define INV_LOW_LINEAR_TRANSFORM_MOVE 6
#define UP_LINEAR_TRANSFORM_MOVE 7
#define INV_UP_LINEAR_TRANSFORM_MOVE 8

#define UPLT 1
#define LOWLT 0
namespace dd {

    typedef struct Move	{
		short index; //执行操作的那一level的
		short pos; //执行操作后的所在的位置
		short optype; //0表示swap，1表示lt，2表示反lt
        unsigned int ddsize; //执行操作后的size
		Move() {
			this->index = -1;
			this->pos = -1;
			this->optype = -1;
			this->ddsize = 0;
		}
		bool operator==(const Move b) const  
		{ 
			return this->index == b.index && this->optype == b.optype && this->ddsize == b.ddsize && this->pos == b.pos;  
		}  
	}Move;

	typedef struct xorNode {
		unsigned short tvar; //上层 dd变量 
		unsigned short cvar; //下层 dd变量
		xorNode(unsigned short hight_DDqubit, unsigned short low_DDqubit) {
			this->tvar = hight_DDqubit;
			this->cvar = low_DDqubit;
		}
	}xorNode;

    enum class dynBuildStrat {
		None, Sifting, LinearSift, SiftingConv, LinearSiftConv, SiftingThenLS, GoodOrderLS, InOrder
	};


}

#endif
