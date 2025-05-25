#define LEFT {1.0f, 0.0f, 0.0f}
#define UP {0.0f, 1.0f, 0.0f}
#define FORWARD {0.0f, 0.0f, 1.0f}

vector extract_forward_from_transform(matrix transform)
{
    return UP * (matrix3)transform;
}