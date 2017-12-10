import os
import subprocess


class Menu(object):
    '''
    combines a menu definition and action dispatcher.  See test example at end of class.
    '''
    keyseparator = ')'
    optionseparator = '  '

    def __init__(self):
        '''
        initialize empty menu object
        usage
            menu = Menu()
        '''
        self.all = {}
        self.query = {}
        self.title = {}
        self.dispatcher = {}
        self.response = ''
        self.menuid = ''

        return None

    def clear(self):
        '''
        clear the windows screen
        :return: True
        '''
        # print('\033[H\033[J')
        # print('\f')
        # os.system('cls')
        print('\n' * 5)
        return True

    def add(self, id, fields):
        '''
        add a menu to the object
        :param id: name of menu
        :param fields: dictionary with the keys and strings for each option
        :return: True
        '''
        newmenu = {}
        query = ''
        for k in fields:
            newmenu[k] = fields[k].replace(k, '{0}{1}'.format(k.upper(), Menu.keyseparator))
        self.all[id] = newmenu

        query = Menu.optionseparator.join(newmenu.values())
        self.query[id] = query
        # print('id:', id, 'query:', query)

        return True

    def ask(self, id, indent):
        '''
        display the menu and get the response
        :param id: name of the menu
        :param indent: number of spaces to indent the query string
        :return: character read from keyboard
        usage
            while True:
                response = menu.ask('top', 2)
        '''
        while True:
            self.clear()
            response = input('\n{0}{1}: '.format(' ' * indent, self.query[id]))
            try:
                # CR only, repeat query
                response = response.upper()[0]
            except IndexError:
                continue

            if response in self.all[id]:
                self.menuid = id
                self.response = response
                return response
            #                    break
            else:
                print('Menu:{0} unknown option {1}\n'.format(id, response))
                continue

    def addTitle(self, id, title):
        '''
        add a descriptinve title to the menu
        :param id: name of menu
        :param title: the printable title for this menu, to be printed above the table
        :return: True
        '''
        self.title[id] = title
        return True

    def addDispatch(self, id, dispatch):

        '''
        Add a dictionary of dispatcher functions for the menu options.
        Must be indexed the same as all and query attributes
        :param id: name of menu
        :param dispatch: dictionary of functions to dispatch when menu item is selected
        :return: True
        '''
        self.dispatcher[id] = dispatch
        return True


    def dispatch(self):
        '''
        call the dispatch function for the current response
        Using the internal value of response and menuide, dispatch the sleected function
        '''
        func = self.dispatcher[self.menuid][self.response]
        # func()
        return func


    def testA(self):
        '''
        dummy function for testing
       '''
        print('function testA')
        return None


    def testQ(self):
        '''
        dummy function for testing
       '''
        print('function testQ')
        return None


if __name__ == '__main__':
    print('Menu::Testing ')

    menu = Menu()
    menu.add('top', {'A': 'Add', 'Q': 'Quit'})
    menu.addDispatch('top', {'A': menu.testA, 'Q': menu.testQ})
    while True:
        response = menu.ask('top', 2)
        print('response is', response)
        menu.dispatch()
        if menu.response == 'Q':
            break

    print('Thank you')
