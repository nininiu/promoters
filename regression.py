"""
View more, visit my tutorial page: https://mofanpy.com/tutorials/
My Youtube Channel: https://www.youtube.com/user/MorvanZhou

Dependencies:
torch: 0.4
matplotlib
"""
import torch
import torch.utils.data as Data
import torch.nn.functional as F
import matplotlib.pyplot as plt
import numpy as np


LR = 0.01
BATCH_SIZE = 32
EPOCH = 1200

f1 = False
# torch.manual_seed(1)    # reproducible
def open_npy_list(input_file):
    data=np.load(input_file,allow_pickle=True)
    data = data.tolist()
    return data

x = open_npy_list('xpssm.npy')
y = open_npy_list('y.npy')
x = torch.FloatTensor(np.array(x))
y = torch.FloatTensor(np.array(y))

# put dateset into torch dataset
torch_dataset = Data.TensorDataset(x, y)
loader = Data.DataLoader(dataset=torch_dataset, batch_size=BATCH_SIZE, shuffle=True, num_workers=2,)

# x = torch.unsqueeze(torch.linspace(-1, 1, 100), dim=1)  # x data (tensor), shape=(100, 1)
# y = x.pow(2) + 0.2*torch.rand(x.size())                 # noisy y data (tensor), shape=(100, 1)

# torch can only train on Variable, so convert them to Variable
# The code below is deprecated in Pytorch 0.4. Now, autograd directly supports tensors
# x, y = Variable(x), Variable(y)

# plt.scatter(x.data.numpy(), y.data.numpy())
# plt.show()

if f1:
    class Net(torch.nn.Module):
        def __init__(self, n_feature, n_hidden, n_output):
            super(Net, self).__init__()
            self.hidden = torch.nn.Linear(n_feature, n_hidden)   # hidden layer
            self.predict = torch.nn.Linear(n_hidden, n_output)   # output layer

        def forward(self, x):
            x = F.sigmoid(self.hidden(x))      # activation function for hidden layer
            x = self.predict(x)             # linear output
            return x
    def original():
        net = Net(n_feature=5, n_hidden=10, n_output=1)     # define the network
        print(net)  # net architecture

    def save():
        net1 = torch.nn.Sequential(
            torch.nn.Linear(5, 10),
            torch.nn.Sigmoid(),
            torch.nn.Linear(10, 1)
        )
        print(net1)  # net architecture

        optimizer = torch.optim.SGD(net1.parameters(), lr=0.2)
        loss_func = torch.nn.MSELoss()  # this is for regression mean squared loss

        for t in range(500000):
            prediction = net1(x)     # input x and predict based on x
            loss = loss_func(prediction, y)     # must be (1. nn output, 2. target)

            optimizer.zero_grad()   # clear gradients for next train
            loss.backward()         # backpropagation, compute gradients
            optimizer.step()        # apply gradients

            if t % 5000 == 0:
                print(loss)
        # plot result
        plt.figure(1, figsize=(10, 3))
        plt.subplot(131)
        plt.title('Net1')
        plt.scatter(y.data.numpy(), prediction.data.numpy())

        # 2 ways to save the net
        torch.save(net1, 'net.pkl')  # save entire net
        torch.save(net1.state_dict(), 'net_params.pkl')   # save only the parameters
        print(net1.state_dict())

    def restore_net():
        # restore entire net1 to net2
        net2 = torch.load('net.pkl')
        prediction = net2(x)
        print(net2)

        # plot result
        plt.subplot(132)
        plt.title('Net2')
        plt.scatter(y.data.numpy(), prediction.data.numpy())


    def restore_params():
        # restore only the parameters in net1 to net3
        net3 = torch.nn.Sequential(
            torch.nn.Linear(5, 10),
            torch.nn.ReLU(),
            torch.nn.Linear(10, 1)
        )
        print(net3)

        # copy net1's parameters into net3
        net3.load_state_dict(torch.load('net_params.pkl'))
        prediction = net3(x)

        # plot result
        plt.subplot(133)
        plt.title('Net3')
        plt.scatter(y.data.numpy(), prediction.data.numpy())
        plt.savefig('./save_reload.jpg')
        plt.show()


    # save net1
    save()

    # restore entire net (may slow)
    restore_net()

    # restore only the net parameters
    restore_params()


# default network
class Net(torch.nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        self.hidden = torch.nn.Linear(5, 20)   # hidden layer
        self.predict = torch.nn.Linear(20, 1)   # output layer

    def forward(self, x):
        x = F.sigmoid(self.hidden(x))      # activation function for hidden layer
        x = self.predict(x)             # linear output
        return x

if __name__ == '__main__':
    # different nets
    net_SGD         = Net()
    net_Momentum    = Net()
    net_RMSprop     = Net()
    net_Adam        = Net()
    nets = [net_SGD, net_Momentum, net_RMSprop, net_Adam]

    # different optimizers
    opt_SGD         = torch.optim.SGD(net_SGD.parameters(), lr=LR)
    opt_Momentum    = torch.optim.SGD(net_Momentum.parameters(), lr=LR, momentum=0.8)
    opt_RMSprop     = torch.optim.RMSprop(net_RMSprop.parameters(), lr=LR, alpha=0.9)
    opt_Adam        = torch.optim.Adam(net_Adam.parameters(), lr=LR, betas=(0.9, 0.99))
    optimizers = [opt_SGD, opt_Momentum, opt_RMSprop, opt_Adam]

    loss_func = torch.nn.MSELoss()
    losses_his = [[], [], [], []]   # record loss

    # training
    for epoch in range(EPOCH):
        print('Epoch: ', epoch)
        for step, (b_x, b_y) in enumerate(loader):          # for each training step
            for net, opt, l_his in zip(nets, optimizers, losses_his):
                output = net(b_x)              # get output for every net
                loss = loss_func(output, b_y)  # compute loss for every net
                opt.zero_grad()                # clear gradients for next train
                loss.backward()                # backpropagation, compute gradients
                opt.step()                     # apply gradients
                l_his.append(loss.data.numpy())     # loss recoder
    for item in losses_his:
        print(item[-1])

    labels = ['SGD', 'Momentum', 'RMSprop', 'Adam']
    for i, l_his in enumerate(losses_his):
        plt.plot(l_his, label=labels[i])
    plt.legend(loc='best')
    plt.xlabel('Steps')
    plt.ylabel('Loss')
    plt.ylim((0, 1))
    plt.savefig('./optimizer.jpg')
    plt.show()
